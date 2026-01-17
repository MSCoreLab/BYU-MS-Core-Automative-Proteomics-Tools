#!/usr/bin/env python3
"""
FASTA File Processor (GUI)

Features:
1. Filter FASTA entries by header text
   - Removes any entries whose HEADER (the line starting with '>') matches user-provided patterns.
   - Patterns can be plain substrings (default) or regular expressions.
   - Options for case sensitivity and saving a removal report.

2. Merge multiple FASTA files
   - Combine multiple FASTA files into a single output file.
   - Optional deduplication based on headers or full sequences.
   - Adds optional prefix to headers from each file.

Example:
  If the pattern is '##', entries like:
    >##gnl|ECOLI|ABC-MONOMER ##L-methionine ...
  will be removed from the output.

Architecture:
  - FastaEntry: Data class representing a single FASTA entry
  - FastaReader: Iterator for reading FASTA files
  - PatternMatcher: Abstract base class for matching strategies
    - SubstringMatcher: Substring-based matching
    - RegexMatcher: Regular expression matching
  - FastaFilter: Handles filtering operations with statistics
  - DeduplicationStrategy: Abstract base class for deduplication
    - NoDuplication: No deduplication
    - HeaderDeduplication: Deduplicate by header
    - SequenceDeduplication: Deduplicate by sequence hash
  - FastaMerger: Handles merge operations with statistics
  - App: Tkinter GUI application
"""

import hashlib
import re
import tkinter as tk
from abc import ABC, abstractmethod
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path
from tkinter import filedialog, messagebox, ttk


@dataclass
class FastaEntry:
    """Represents a single FASTA entry with header and sequence.

    Attributes:
        header: The header line including leading '>'
        sequence_lines: List of sequence lines (preserves original formatting)
    """
    header: str
    sequence_lines: list[str]

    @property
    def header_text(self) -> str:
        """Return header without leading '>'.

        Returns:
            Header text with '>' prefix removed
        """
        return self.header.lstrip(">")

    @property
    def sequence(self) -> str:
        """Return joined sequence.

        Returns:
            Complete sequence as single string
        """
        return "".join(self.sequence_lines)

    @property
    def sequence_hash(self) -> str:
        """Return MD5 hash of sequence for deduplication.

        Returns:
            Hexadecimal MD5 hash string

        Note:
            usedforsecurity=False since this is for data deduplication,
            not cryptographic purposes
        """
        return hashlib.md5(self.sequence.encode(), usedforsecurity=False).hexdigest()

    def write_to_file(self, file_handle, prefix: str = ""):
        """Write entry to an open file handle with optional prefix.

        Args:
            file_handle: Open file object for writing
            prefix: Optional prefix to add to header (e.g., '[filename]')
        """
        # * Write header with optional prefix
        if prefix:
            file_handle.write(f">{prefix}{self.header_text}\n")
        else:
            file_handle.write(f"{self.header}\n")

        # * Write sequence lines preserving original formatting
        for line in self.sequence_lines:
            file_handle.write(f"{line}\n")


class FastaReader:
    """Reads and iterates over FASTA entries from a file.

    Uses iterator pattern for memory-efficient processing of large FASTA files.
    Handles malformed files gracefully by replacing invalid UTF-8 characters.
    """

    def __init__(self, file_path: Path):
        """Initialize reader with file path.

        Args:
            file_path: Path object pointing to FASTA file
        """
        self.file_path = file_path

    def __iter__(self) -> Iterator[FastaEntry]:
        """Yield FastaEntry objects from the file.

        Yields:
            FastaEntry: Individual FASTA entries (header + sequence lines)

        Note:
            - Lines starting with '>' are treated as headers
            - Sequence lines are accumulated until next header
            - Final entry is yielded after file ends
        """
        with self.file_path.open("r", encoding="utf-8", errors="replace") as f:
            header = None
            seq_lines = []

            for line in f:
                line = line.rstrip("\n")

                # ! New header found - yield previous entry if exists
                if line.startswith(">"):
                    if header is not None:
                        yield FastaEntry(header, seq_lines)
                    header = line
                    seq_lines = []
                # * Sequence line - append to current entry
                elif header is not None:
                    seq_lines.append(line)

            # ! Yield final entry if file didn't end with newline
            if header is not None:
                yield FastaEntry(header, seq_lines)


class PatternMatcher(ABC):
    """Abstract base class for pattern matching strategies."""

    @abstractmethod
    def matches(self, text: str) -> bool:
        """Check if text matches any pattern."""
        pass


class SubstringMatcher(PatternMatcher):
    """Matches patterns as substrings.

    Fast substring matching using Python's built-in 'in' operator.
    Normalizes both patterns and text when case-insensitive.
    """

    def __init__(self, patterns: list[str], case_sensitive: bool = False):
        """Initialize matcher with patterns.

        Args:
            patterns: List of substring patterns to match
            case_sensitive: Whether to perform case-sensitive matching
        """
        # * Normalize patterns once at initialization for performance
        self.patterns = patterns if case_sensitive else [p.lower() for p in patterns]
        self.case_sensitive = case_sensitive

    def matches(self, text: str) -> bool:
        """Check if text contains any pattern.

        Args:
            text: Text to search for patterns

        Returns:
            True if any pattern is found in text
        """
        normalized_text = text if self.case_sensitive else text.lower()
        return any(pattern in normalized_text for pattern in self.patterns)


class RegexMatcher(PatternMatcher):
    """Matches patterns as regular expressions.

    Compiles regex patterns at initialization for better performance.
    Validates regex syntax and provides helpful error messages.
    """

    def __init__(self, patterns: list[str], case_sensitive: bool = False):
        """Initialize matcher with regex patterns.

        Args:
            patterns: List of regex pattern strings
            case_sensitive: Whether to perform case-sensitive matching

        Raises:
            ValueError: If any pattern has invalid regex syntax
        """
        flags = 0 if case_sensitive else re.IGNORECASE
        try:
            # ! Pre-compile all regex patterns for performance
            self.regexes = [re.compile(p, flags) for p in patterns]
        except re.error as e:
            # ! Provide user-friendly error message for invalid regex
            raise ValueError(f"Invalid regular expression: {e}") from e

    def matches(self, text: str) -> bool:
        """Check if text matches any regex pattern.

        Args:
            text: Text to search for pattern matches

        Returns:
            True if any pattern matches the text
        """
        return any(regex.search(text) for regex in self.regexes)


@dataclass
class FilterStats:
    """Statistics from filtering operation."""
    kept: int = 0
    removed: int = 0
    removed_headers: list[str] = field(default_factory=list)


@dataclass
class MergeStats:
    """Statistics from merge operation."""
    total_entries: int = 0
    written_entries: int = 0
    skipped_duplicates: int = 0
    file_stats: dict[str, tuple[int, int]] = field(default_factory=dict)


class FastaFilter:
    """Handles FASTA filtering operations.

    Removes entries whose headers match specified patterns using either
    substring or regex matching. Uses Strategy pattern for flexible matching.
    """

    def __init__(self, patterns: list[str], use_regex: bool = False, case_sensitive: bool = False):
        """Initialize filter with patterns and matching options.

        Args:
            patterns: List of patterns to match against headers
            use_regex: Whether to treat patterns as regular expressions
            case_sensitive: Whether matching should be case-sensitive

        Raises:
            ValueError: If patterns list is empty
        """
        if not patterns:
            raise ValueError("Please provide at least one pattern to match.")

        # * Use Strategy pattern to select matcher implementation
        self.matcher = (
            RegexMatcher(patterns, case_sensitive) if use_regex
            else SubstringMatcher(patterns, case_sensitive)
        )
        # * Store configuration for report generation
        self.patterns = patterns
        self.use_regex = use_regex
        self.case_sensitive = case_sensitive

    def filter_file(self, input_path: Path, output_path: Path) -> FilterStats:
        """Filter FASTA file, removing entries matching patterns.

        Args:
            input_path: Path to input FASTA file
            output_path: Path for filtered output file

        Returns:
            FilterStats object containing operation statistics

        Note:
            - Only headers are checked; sequences pass through unchanged
            - Memory efficient: processes one entry at a time
        """
        stats = FilterStats()
        reader = FastaReader(input_path)

        with output_path.open("w", encoding="utf-8") as fout:
            for entry in reader:
                # ! Check if header matches any removal pattern
                if self.matcher.matches(entry.header_text):
                    stats.removed += 1
                    stats.removed_headers.append(entry.header)
                # * Keep entry - write to output file
                else:
                    stats.kept += 1
                    entry.write_to_file(fout)

        return stats

    def save_report(self, input_path: Path, output_path: Path, stats: FilterStats, report_path: Path):
        """Save filtering report to file."""
        with report_path.open("w", encoding="utf-8") as rep:
            rep.write(f"Input:  {input_path}\n")
            rep.write(f"Output: {output_path}\n")
            rep.write(f"Patterns: {self.patterns}\n")
            rep.write(f"Regex: {self.use_regex}\n")
            rep.write(f"Case sensitive: {self.case_sensitive}\n\n")
            rep.write(f"Kept entries:    {stats.kept}\n")
            rep.write(f"Removed entries: {stats.removed}\n\n")
            if stats.removed_headers:
                rep.write("Removed headers:\n")
                for header in stats.removed_headers:
                    rep.write(f"{header}\n")


class DeduplicationStrategy(ABC):
    """Abstract base class for deduplication strategies."""

    @abstractmethod
    def is_duplicate(self, entry: FastaEntry) -> bool:
        """Check if entry is a duplicate."""
        pass


class NoDuplication(DeduplicationStrategy):
    """No deduplication - keep all entries.

    Null Object pattern: always returns False for is_duplicate.
    Used when user doesn't want any deduplication.
    """

    def is_duplicate(self, entry: FastaEntry) -> bool:
        """Always returns False (never considers entries as duplicates).

        Args:
            _entry: Unused parameter (prefixed with _ to indicate intentionally unused)

        Returns:
            Always False
        """
        _ = entry
        return False


class HeaderDeduplication(DeduplicationStrategy):
    """Deduplicate by header text.

    Maintains a set of seen headers. First occurrence is kept,
    subsequent entries with same header are marked as duplicates.
    """

    def __init__(self):
        """Initialize with empty set of seen headers."""
        self.seen_headers: set[str] = set()

    def is_duplicate(self, entry: FastaEntry) -> bool:
        """Check if entry's header was seen before.

        Args:
            entry: FASTA entry to check

        Returns:
            True if header was seen before, False otherwise

        Side Effect:
            Adds header to seen_headers set on first occurrence
        """
        header_text = entry.header_text
        # * Check if header already seen
        if header_text in self.seen_headers:
            return True
        # * First occurrence - add to set and keep
        self.seen_headers.add(header_text)
        return False


class SequenceDeduplication(DeduplicationStrategy):
    """Deduplicate by sequence hash.

    Uses MD5 hashing for efficient sequence comparison.
    First occurrence is kept, subsequent entries with identical
    sequences are marked as duplicates.

    Note:
        More memory-efficient than storing full sequences.
        Hash collisions are extremely unlikely for typical FASTA files.
    """

    def __init__(self):
        """Initialize with empty set of seen sequence hashes."""
        self.seen_hashes: set[str] = set()

    def is_duplicate(self, entry: FastaEntry) -> bool:
        """Check if entry's sequence was seen before.

        Args:
            entry: FASTA entry to check

        Returns:
            True if sequence hash was seen before, False otherwise

        Side Effect:
            Adds sequence hash to seen_hashes set on first occurrence
        """
        seq_hash = entry.sequence_hash
        # * Check if sequence hash already seen
        if seq_hash in self.seen_hashes:
            return True
        # * First occurrence - add hash to set and keep
        self.seen_hashes.add(seq_hash)
        return False


class FastaMerger:
    """Handles FASTA merging operations."""

    DEDUP_STRATEGIES = {
        "none": NoDuplication,
        "header": HeaderDeduplication,
        "sequence": SequenceDeduplication,
    }

    def __init__(self, deduplicate: str = "none", add_prefix: bool = False):
        if deduplicate not in self.DEDUP_STRATEGIES:
            raise ValueError(f"Invalid deduplication mode: {deduplicate}")

        self.dedup_strategy = self.DEDUP_STRATEGIES[deduplicate]()
        self.add_prefix = add_prefix
        self.deduplicate_mode = deduplicate

    def merge_files(self, input_paths: list[Path], output_path: Path) -> MergeStats:
        """Merge multiple FASTA files into one."""
        if not input_paths:
            raise ValueError("Please provide at least one input file.")

        stats = MergeStats()

        with output_path.open("w", encoding="utf-8") as fout:
            for input_path in input_paths:
                if not input_path.exists():
                    raise FileNotFoundError(f"Input file not found: {input_path}")

                file_stats = self._process_file(input_path, fout, stats)
                stats.file_stats[input_path.name] = file_stats

        return stats

    def _process_file(self, input_path: Path, output_handle, stats: MergeStats) -> tuple[int, int]:
        """Process a single file during merge.

        Args:
            input_path: Path to input FASTA file
            output_handle: Open file handle for merged output
            stats: MergeStats object to update

        Returns:
            Tuple of (total_entries_in_file, entries_written_from_file)
        """
        file_total = 0
        file_written = 0
        # * Generate prefix if enabled (e.g., '[filename]')
        prefix = f"[{input_path.stem}]" if self.add_prefix else ""

        reader = FastaReader(input_path)
        for entry in reader:
            file_total += 1
            stats.total_entries += 1

            # * Check for duplicates using selected strategy
            if self.dedup_strategy.is_duplicate(entry):
                stats.skipped_duplicates += 1
            else:
                # * Write entry with optional prefix
                entry.write_to_file(output_handle, prefix)
                stats.written_entries += 1
                file_written += 1

        return file_total, file_written

    def save_report(self, output_path: Path, stats: MergeStats, report_path: Path):
        """Save merge report to file."""
        with report_path.open("w", encoding="utf-8") as rep:
            rep.write("FASTA Merge Report\n")
            rep.write("=" * 50 + "\n\n")
            rep.write(f"Output file: {output_path}\n")
            rep.write(f"Deduplication mode: {self.deduplicate_mode}\n")
            rep.write(f"Add file prefix: {self.add_prefix}\n\n")
            rep.write(f"Total entries processed: {stats.total_entries}\n")
            rep.write(f"Total entries written: {stats.written_entries}\n")
            rep.write(f"Duplicate entries skipped: {stats.skipped_duplicates}\n\n")
            rep.write("Per-file statistics:\n")
            rep.write("-" * 50 + "\n")
            for filename, (total, written) in stats.file_stats.items():
                rep.write(f"{filename}: {total} entries, {written} written\n")


class App(tk.Tk):
    # Class constants
    FASTA_FILETYPES = [
        ("FASTA files", "*.fasta *.fa *.faa *.fna *.fas *.fsa"),
        ("All files", "*.*"),
    ]
    PADDING = {"padx": 10, "pady": 8}

    # Dark mode colors
    DARK_BG = "#1e1e1e"
    DARK_FG = "#e0e0e0"
    DARK_ACCENT = "#3c3c3c"
    DARK_HIGHLIGHT = "#007acc"

    def __init__(self):
        """Initialize the FASTA File Processor GUI application."""
        super().__init__()
        self.title("FASTA File Processor")
        self.geometry("800x500")
        self.minsize(700, 450)

        # * Apply dark mode theme
        self._setup_dark_theme()

        # * Create notebook for tabs
        self.notebook = ttk.Notebook(self)
        self.notebook.pack(fill="both", expand=True, padx=10, pady=10)

        # * Create filter tab
        self.filter_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.filter_frame, text="Filter FASTA")

        # * Create merge tab
        self.merge_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.merge_frame, text="Merge FASTA Files")

        # * Initialize filter tab variables
        self.input_var = tk.StringVar()
        self.output_var = tk.StringVar()
        self.patterns_var = tk.StringVar(value="##")  # Default pattern
        self.regex_var = tk.BooleanVar(value=False)
        self.case_var = tk.BooleanVar(value=False)
        self.report_var = tk.BooleanVar(value=True)  # Generate reports by default

        # * Initialize merge tab variables
        self.merge_files: list[Path] = []  # List of Path objects to merge
        self.merge_output_var = tk.StringVar()
        self.dedupe_var = tk.StringVar(value="none")  # No deduplication by default
        self.prefix_var = tk.BooleanVar(value=False)
        self.merge_report_var = tk.BooleanVar(value=True)

        # * Build UI for both tabs
        self._build_filter_ui()
        self._build_merge_ui()

    def _setup_dark_theme(self):
        """Configure dark mode styling for ttk widgets."""
        self.root = self
        self.root.configure(bg=self.DARK_BG)

        style = ttk.Style()
        style.theme_use("clam")

        # Configure dark colors for all ttk widgets
        style.configure(".", background=self.DARK_BG, foreground=self.DARK_FG)
        style.configure("TFrame", background=self.DARK_BG)
        style.configure("TLabel", background=self.DARK_BG, foreground=self.DARK_FG)
        style.configure(
            "TButton",
            background=self.DARK_ACCENT,
            foreground=self.DARK_FG,
            borderwidth=1,
            focuscolor=self.DARK_HIGHLIGHT,
        )
        style.map("TButton", background=[("active", self.DARK_HIGHLIGHT)])
        style.configure("TNotebook", background=self.DARK_BG, borderwidth=0)
        style.configure(
            "TNotebook.Tab", background=self.DARK_ACCENT, foreground=self.DARK_FG, padding=[10, 5]
        )
        style.map(
            "TNotebook.Tab",
            background=[("selected", self.DARK_HIGHLIGHT)],
            foreground=[("selected", "white")],
        )
        style.configure("TCheckbutton", background=self.DARK_BG, foreground=self.DARK_FG)
        style.configure("TRadiobutton", background=self.DARK_BG, foreground=self.DARK_FG)

    def _build_filter_ui(self):
        frm = ttk.Frame(self.filter_frame)
        frm.pack(fill="both", expand=True)

        # Row 1: Input
        row1 = ttk.Frame(frm)
        row1.pack(fill="x")
        ttk.Label(row1, text="Input FASTA:").pack(side="left")
        input_entry = tk.Entry(
            row1,
            textvariable=self.input_var,
            bg=self.DARK_ACCENT,
            fg=self.DARK_FG,
            insertbackground=self.DARK_FG,
        )
        input_entry.pack(side="left", fill="x", expand=True, padx=6)
        ttk.Button(row1, text="Browse…", command=self.choose_input).pack(side="left")

        # Row 2: Output
        row2 = ttk.Frame(frm)
        row2.pack(fill="x")
        ttk.Label(row2, text="Output FASTA:").pack(side="left")
        output_entry = tk.Entry(
            row2,
            textvariable=self.output_var,
            bg=self.DARK_ACCENT,
            fg=self.DARK_FG,
            insertbackground=self.DARK_FG,
        )
        output_entry.pack(side="left", fill="x", expand=True, padx=6)
        ttk.Button(row2, text="Browse…", command=self.choose_output).pack(side="left")

        # Row 3: Patterns
        row3 = ttk.Frame(frm)
        row3.pack(fill="x")
        ttk.Label(row3, text="Patterns (comma-separated):").pack(side="left")
        pattern_entry = tk.Entry(
            row3,
            textvariable=self.patterns_var,
            bg=self.DARK_ACCENT,
            fg=self.DARK_FG,
            insertbackground=self.DARK_FG,
        )
        pattern_entry.pack(side="left", fill="x", expand=True, padx=6)

        # Row 4: Options
        row4 = ttk.Frame(frm)
        row4.pack(fill="x")
        ttk.Checkbutton(row4, text="Use regular expressions", variable=self.regex_var).pack(
            side="left"
        )
        ttk.Checkbutton(row4, text="Case sensitive", variable=self.case_var).pack(side="left")
        ttk.Checkbutton(row4, text="Save removal report (.txt)", variable=self.report_var).pack(
            side="left"
        )

        # Row 5: Actions
        row5 = ttk.Frame(frm)
        row5.pack(fill="x")
        ttk.Button(row5, text="Run Filter", command=self.run_filter).pack(side="left")

        # Help text
        help_txt = (
            "Tips:\n"
            "• Enter one or more patterns separated by commas. Example: ##,gnl|ECOLI|ABC\n"
            "• If 'Use regular expressions' is ON, each pattern is a regex.\n"
            "• Only headers starting with '>' are checked; sequences are preserved as-is for kept entries.\n"
        )
        row6 = ttk.Frame(frm)
        row6.pack(fill="both", expand=True)
        ttk.Label(row6, text=help_txt, justify="left").pack(anchor="w")

    def _build_merge_ui(self):
        frm = ttk.Frame(self.merge_frame)
        frm.pack(fill="both", expand=True)

        # Row 1: File list with scrollbar
        row1 = ttk.Frame(frm)
        row1.pack(fill="both", expand=True)
        ttk.Label(row1, text="Input FASTA files to merge:").pack(anchor="w")

        list_frame = ttk.Frame(row1)
        list_frame.pack(fill="both", expand=True, pady=4)

        scrollbar = ttk.Scrollbar(list_frame)
        scrollbar.pack(side="right", fill="y")

        self.merge_listbox = tk.Listbox(
            list_frame,
            yscrollcommand=scrollbar.set,
            height=6,
            bg=self.DARK_ACCENT,
            fg=self.DARK_FG,
            selectbackground=self.DARK_HIGHLIGHT,
            selectforeground="white",
            borderwidth=0,
            highlightthickness=1,
            highlightcolor=self.DARK_HIGHLIGHT,
        )
        self.merge_listbox.pack(side="left", fill="both", expand=True)
        scrollbar.config(command=self.merge_listbox.yview)

        # Buttons for file list management
        btn_frame = ttk.Frame(row1)
        btn_frame.pack(fill="x", pady=4)
        ttk.Button(btn_frame, text="Add Files…", command=self.add_merge_files).pack(
            side="left", padx=2
        )
        ttk.Button(btn_frame, text="Remove Selected", command=self.remove_merge_file).pack(
            side="left", padx=2
        )
        ttk.Button(btn_frame, text="Clear All", command=self.clear_merge_files).pack(
            side="left", padx=2
        )

        # Row 2: Output
        row2 = ttk.Frame(frm)
        row2.pack(fill="x")
        ttk.Label(row2, text="Output merged FASTA:").pack(side="left")
        merge_output_entry = tk.Entry(
            row2,
            textvariable=self.merge_output_var,
            bg=self.DARK_ACCENT,
            fg=self.DARK_FG,
            insertbackground=self.DARK_FG,
        )
        merge_output_entry.pack(side="left", fill="x", expand=True, padx=6)
        ttk.Button(row2, text="Browse…", command=self.choose_merge_output).pack(side="left")

        # Row 3: Options
        row3 = ttk.Frame(frm)
        row3.pack(fill="x")

        ttk.Label(row3, text="Deduplication:").pack(side="left", padx=(0, 6))
        ttk.Radiobutton(row3, text="None", variable=self.dedupe_var, value="none").pack(side="left")
        ttk.Radiobutton(row3, text="By Header", variable=self.dedupe_var, value="header").pack(
            side="left"
        )
        ttk.Radiobutton(row3, text="By Sequence", variable=self.dedupe_var, value="sequence").pack(
            side="left"
        )

        # Row 4: More options
        row4 = ttk.Frame(frm)
        row4.pack(fill="x")
        ttk.Checkbutton(row4, text="Add file prefix to headers", variable=self.prefix_var).pack(
            side="left"
        )
        ttk.Checkbutton(row4, text="Save merge report (.txt)", variable=self.merge_report_var).pack(
            side="left"
        )

        # Row 5: Actions
        row5 = ttk.Frame(frm)
        row5.pack(fill="x")
        ttk.Button(row5, text="Run Merge", command=self.run_merge).pack(side="left")

        # Help text
        help_txt = (
            "Tips:\n"
            "• Add multiple FASTA files to merge them into one combined file.\n"
            "• Deduplication: 'None' keeps all entries, 'By Header' removes duplicate headers,\n"
            "  'By Sequence' removes duplicate sequences.\n"
            "• File prefix adds source filename to each header (e.g., >[file1]original_header).\n"
        )
        row6 = ttk.Frame(frm)
        row6.pack(fill="x")
        ttk.Label(row6, text=help_txt, justify="left").pack(anchor="w")

    def choose_input(self):
        """Open file dialog to select input FASTA file.

        Side Effects:
            - Sets input_var to selected file path
            - Auto-suggests output filename if not already set
        """
        path = filedialog.askopenfilename(
            title="Choose FASTA file",
            filetypes=self.FASTA_FILETYPES,
        )
        if path:
            self.input_var.set(path)
            # * Suggest an output filename next to input
            in_path = Path(path)
            suggested = in_path.with_suffix(in_path.suffix + ".filtered.fasta")
            # * Only auto-suggest if output not already set
            if not self.output_var.get():
                self.output_var.set(str(suggested))

    def choose_output(self):
        path = filedialog.asksaveasfilename(
            title="Save filtered FASTA as…",
            defaultextension=".fasta",
            filetypes=self.FASTA_FILETYPES,
        )
        if path:
            self.output_var.set(path)

    def add_merge_files(self):
        """Open file dialog to select FASTA files for merging.

        Allows multiple file selection. Prevents duplicate entries.
        Auto-suggests output filename if not already set.

        Side Effects:
            - Adds selected files to merge_files list
            - Updates merge_listbox display
            - Auto-suggests merge output path on first addition
        """
        paths = filedialog.askopenfilenames(
            title="Choose FASTA files to merge",
            filetypes=self.FASTA_FILETYPES,
        )
        if paths:
            for path in paths:
                p = Path(path)
                # * Prevent duplicate file additions
                if p not in self.merge_files:
                    self.merge_files.append(p)
                    self.merge_listbox.insert(tk.END, p.name)

            # * Suggest output if not set (use parent dir of first file)
            if not self.merge_output_var.get() and self.merge_files:
                first = self.merge_files[0]
                suggested = first.parent / "merged_output.fasta"
                self.merge_output_var.set(str(suggested))

    def remove_merge_file(self):
        selection = self.merge_listbox.curselection()
        if selection:
            idx = selection[0]
            self.merge_listbox.delete(idx)
            del self.merge_files[idx]

    def clear_merge_files(self):
        self.merge_listbox.delete(0, tk.END)
        self.merge_files.clear()

    def choose_merge_output(self):
        path = filedialog.asksaveasfilename(
            title="Save merged FASTA as…",
            defaultextension=".fasta",
            filetypes=self.FASTA_FILETYPES,
        )
        if path:
            self.merge_output_var.set(path)

    def _parse_patterns(self, patterns_str: str) -> list[str]:
        """Parse comma-separated patterns and validate."""
        patterns = [p.strip() for p in patterns_str.split(",") if p.strip()]
        if not patterns:
            raise ValueError("Please enter at least one pattern.")
        return patterns

    def run_filter(self):
        """Execute FASTA filtering operation.

        Validates inputs, creates FastaFilter instance, processes file,
        and displays results. Optionally generates removal report.

        Shows:
            - Success dialog with statistics
            - Error dialog if validation fails or exception occurs
        """
        try:
            # * Validate and expand file paths
            input_path = Path(self.input_var.get()).expanduser()
            output_path = Path(self.output_var.get()).expanduser()

            # ! Validate input file exists
            if not input_path or not input_path.exists():
                messagebox.showerror("Error", "Please choose a valid input FASTA file.")
                return
            # ! Validate output path is set
            if not output_path:
                messagebox.showerror("Error", "Please choose an output FASTA file path.")
                return

            # * Parse and validate patterns
            patterns = self._parse_patterns(self.patterns_var.get().strip())

            # * Create filter and process file
            fasta_filter = FastaFilter(
                patterns=patterns,
                use_regex=self.regex_var.get(),
                case_sensitive=self.case_var.get()
            )

            stats = fasta_filter.filter_file(input_path, output_path)

            # * Save report if requested
            if self.report_var.get():
                report_path = output_path.with_suffix(output_path.suffix + ".removed.txt")
                fasta_filter.save_report(input_path, output_path, stats, report_path)
                msg = f"Done!\nKept entries: {stats.kept}\nRemoved entries: {stats.removed}"
                msg += f"\n\nReport saved to:\n{report_path}"
            else:
                msg = f"Done!\nKept entries: {stats.kept}\nRemoved entries: {stats.removed}"

            messagebox.showinfo("FASTA Header Filter", msg)
        except Exception as e:
            # ! Display user-friendly error message
            messagebox.showerror("Error", str(e))

    def run_merge(self):
        """Execute FASTA merge operation.

        Validates inputs, creates FastaMerger instance, processes files,
        and displays results. Optionally generates merge report.

        Shows:
            - Success dialog with statistics
            - Error dialog if validation fails or exception occurs
        """
        try:
            # ! Validate at least one file selected
            if not self.merge_files:
                messagebox.showerror("Error", "Please add at least one FASTA file to merge.")
                return

            # * Validate and expand output path
            output_path = Path(self.merge_output_var.get()).expanduser()
            if not output_path:
                messagebox.showerror("Error", "Please choose an output FASTA file path.")
                return

            # * Create merger and process files
            merger = FastaMerger(
                deduplicate=self.dedupe_var.get(),
                add_prefix=self.prefix_var.get()
            )

            stats = merger.merge_files(self.merge_files, output_path)

            # * Save report if requested
            if self.merge_report_var.get():
                report_path = output_path.with_suffix(output_path.suffix + ".merge_report.txt")
                merger.save_report(output_path, stats, report_path)
                msg = f"Done!\nTotal entries written: {stats.written_entries}"
                if stats.skipped_duplicates > 0:
                    msg += f"\nDuplicate entries skipped: {stats.skipped_duplicates}"
                msg += f"\n\nReport saved to:\n{report_path}"
            else:
                msg = f"Done!\nTotal entries written: {stats.written_entries}"
                if stats.skipped_duplicates > 0:
                    msg += f"\nDuplicate entries skipped: {stats.skipped_duplicates}"

            messagebox.showinfo("FASTA Merge", msg)
        except Exception as e:
            # ! Display user-friendly error message
            messagebox.showerror("Error", str(e))


def main():
    """Main entry point for FASTA File Processor.

    Supports both GUI and command-line (headless) operation:
    - GUI mode: No arguments provided
    - CLI mode: input, output, and patterns arguments provided

    CLI Usage:
        python filter_fasta_gui.py input.fasta output.fasta "pat1,pat2" [--regex] [--case] [--report]
    """
    import argparse

    # * Setup argument parser for CLI mode
    parser = argparse.ArgumentParser(description="Filter FASTA entries by header patterns.")
    parser.add_argument("input", nargs="?", help="Input FASTA")
    parser.add_argument("output", nargs="?", help="Output FASTA")
    parser.add_argument("patterns", nargs="?", help="Comma-separated patterns")
    parser.add_argument(
        "--regex", action="store_true", help="Treat patterns as regular expressions"
    )
    parser.add_argument("--case", action="store_true", help="Case-sensitive matching")
    parser.add_argument(
        "--report", action="store_true", help="Save a removal report next to output"
    )

    args = parser.parse_args()

    # * Check if running in CLI mode (all required args provided)
    if args.input and args.output and args.patterns:
        # ! CLI mode - process file without GUI
        patterns = [p.strip() for p in args.patterns.split(",") if p.strip()]
        input_path = Path(args.input)
        output_path = Path(args.output)

        # * Create filter and process
        fasta_filter = FastaFilter(
            patterns=patterns,
            use_regex=args.regex,
            case_sensitive=args.case
        )

        stats = fasta_filter.filter_file(input_path, output_path)

        # * Save report if requested
        if args.report:
            report_path = output_path.with_suffix(output_path.suffix + ".removed.txt")
            fasta_filter.save_report(input_path, output_path, stats, report_path)

        # * Print results to stdout
        print(f"Kept entries: {stats.kept}")
        print(f"Removed entries: {stats.removed}")
    else:
        # * GUI mode - launch Tkinter application
        app = App()
        app.mainloop()


if __name__ == "__main__":
    main()
