 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
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
"""

import re
from pathlib import Path
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from typing import Optional

def filter_fasta(
    input_path: Path,
    output_path: Path,
    patterns,
    use_regex: bool = False,
    case_sensitive: bool = False,
    save_report_path: Optional[Path] = None,
):
    """
    Read a FASTA, writing only entries whose headers DO NOT match any of the patterns.

    Parameters
    ----------
    input_path : Path
        Source FASTA file.
    output_path : Path
        Destination FASTA file (will be overwritten).
    patterns : list[str]
        List of patterns. If use_regex=False, they are treated as substrings.
        If use_regex=True, they are compiled as regular expressions.
    use_regex : bool
        If True, treat patterns as regex.
    case_sensitive : bool
        If False, perform case-insensitive matching.
    save_report_path : Path | None
        If provided, write a report listing removed headers.
    """
    if not patterns:
        raise ValueError("Please provide at least one pattern to match.")

    # Prepare matchers
    flags = 0 if case_sensitive else re.IGNORECASE
    if use_regex:
        try:
            regexes = [re.compile(p, flags) for p in patterns]
        except re.error as e:
            raise ValueError(f"Invalid regular expression: {e}") from e

        def matches(header: str) -> bool:
            return any(rx.search(header) for rx in regexes)
    else:
        # Normalize case if needed
        pats = patterns if case_sensitive else [p.lower() for p in patterns]

        def matches(header: str) -> bool:
            h = header if case_sensitive else header.lower()
            return any(p in h for p in pats)

    kept = 0
    removed = 0
    removed_headers = []

    with input_path.open('r', encoding='utf-8', errors='replace') as fin, \
         output_path.open('w', encoding='utf-8') as fout:
        current_entry = []  # list[str] of original lines

        def flush_entry():
            nonlocal kept, removed
            if not current_entry:
                return
            header_line = current_entry[0].rstrip('\n')
            header_text = header_line[1:]  # drop leading '>'
            if matches(header_text):
                removed += 1
                removed_headers.append(header_line)
            else:
                kept += 1
                # Write back exactly as read
                for line in current_entry:
                    fout.write(line if line.endswith('\n') else line + '\n')

        for line in fin:
            if line.startswith('>'):
                # New header -> flush previous entry
                flush_entry()
                current_entry = [line.rstrip('\n')]
            else:
                if current_entry:
                    current_entry.append(line.rstrip('\n'))
                else:
                    # Handle files that start with sequence without header (non-standard)
                    # We'll skip such lines but you could also choose to pass them through.
                    continue
        # Flush last entry
        flush_entry()

    if save_report_path is not None:
        with save_report_path.open('w', encoding='utf-8') as rep:
            rep.write(f"Input:  {input_path}\n")
            rep.write(f"Output: {output_path}\n")
            rep.write(f"Patterns: {patterns}\n")
            rep.write(f"Regex: {use_regex}\n")
            rep.write(f"Case sensitive: {case_sensitive}\n\n")
            rep.write(f"Kept entries:    {kept}\n")
            rep.write(f"Removed entries: {removed}\n\n")
            if removed_headers:
                rep.write("Removed headers:\n")
                for h in removed_headers:
                    rep.write(h + "\n")

    return kept, removed


def merge_fasta_files(
    input_paths: list[Path],
    output_path: Path,
    deduplicate: str = "none",
    add_prefix: bool = False,
    save_report_path: Optional[Path] = None,
):
    """
    Merge multiple FASTA files into one.

    Parameters
    ----------
    input_paths : list[Path]
        List of source FASTA files.
    output_path : Path
        Destination FASTA file (will be overwritten).
    deduplicate : str
        Deduplication mode:
        - "none": Keep all entries from all files
        - "header": Remove entries with duplicate headers (first occurrence kept)
        - "sequence": Remove entries with duplicate sequences (first occurrence kept)
    add_prefix : bool
        If True, add a prefix to headers indicating source file (e.g., "[file1]").
    save_report_path : Path | None
        If provided, write a report with merge statistics.
    """
    if not input_paths:
        raise ValueError("Please provide at least one input file.")

    seen_headers = set()
    seen_sequences = set()
    total_entries = 0
    written_entries = 0
    skipped_duplicates = 0
    file_stats = {}  # file -> (total, written)

    with output_path.open('w', encoding='utf-8') as fout:
        for input_path in input_paths:
            if not input_path.exists():
                raise FileNotFoundError(f"Input file not found: {input_path}")
            
            file_total = 0
            file_written = 0
            
            # Get a short prefix from filename (without extension)
            prefix = f"[{input_path.stem}]" if add_prefix else ""
            
            with input_path.open('r', encoding='utf-8', errors='replace') as fin:
                current_header = None
                current_sequence_lines = []
                
                def flush_entry():
                    nonlocal written_entries, skipped_duplicates, file_written
                    if current_header is None:
                        return
                    
                    # Check for duplicates
                    should_write = True
                    
                    if deduplicate == "header":
                        # Extract the actual header text (without '>')
                        header_text = current_header.lstrip('>')
                        if header_text in seen_headers:
                            should_write = False
                            skipped_duplicates += 1
                        else:
                            seen_headers.add(header_text)
                    
                    elif deduplicate == "sequence":
                        # Combine all sequence lines into one string
                        seq = ''.join(current_sequence_lines)
                        if seq in seen_sequences:
                            should_write = False
                            skipped_duplicates += 1
                        else:
                            seen_sequences.add(seq)
                    
                    if should_write:
                        # Write header with optional prefix
                        if prefix and not current_header.startswith('>'):
                            fout.write(f">{prefix}{current_header}\n")
                        elif prefix:
                            fout.write(f">{prefix}{current_header.lstrip('>')}\n")
                        else:
                            fout.write(f"{current_header}\n")
                        
                        # Write sequence lines
                        for seq_line in current_sequence_lines:
                            fout.write(f"{seq_line}\n")
                        
                        written_entries += 1
                        file_written += 1
                
                for line in fin:
                    line = line.rstrip('\n')
                    if line.startswith('>'):
                        # Flush previous entry
                        flush_entry()
                        # Start new entry
                        current_header = line
                        current_sequence_lines = []
                        file_total += 1
                        total_entries += 1
                    else:
                        if current_header is not None:
                            current_sequence_lines.append(line)
                
                # Flush last entry from this file
                flush_entry()
            
            file_stats[input_path.name] = (file_total, file_written)
    
    # Generate report if requested
    if save_report_path is not None:
        with save_report_path.open('w', encoding='utf-8') as rep:
            rep.write("FASTA Merge Report\n")
            rep.write("=" * 50 + "\n\n")
            rep.write(f"Output file: {output_path}\n")
            rep.write(f"Deduplication mode: {deduplicate}\n")
            rep.write(f"Add file prefix: {add_prefix}\n\n")
            rep.write(f"Total entries processed: {total_entries}\n")
            rep.write(f"Total entries written: {written_entries}\n")
            rep.write(f"Duplicate entries skipped: {skipped_duplicates}\n\n")
            rep.write("Per-file statistics:\n")
            rep.write("-" * 50 + "\n")
            for filename, (total, written) in file_stats.items():
                rep.write(f"{filename}: {total} entries, {written} written\n")
    
    return written_entries, skipped_duplicates


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("FASTA File Processor")
        self.geometry("800x500")
        self.minsize(700, 450)

        # Create notebook for tabs
        self.notebook = ttk.Notebook(self)
        self.notebook.pack(fill="both", expand=True, padx=10, pady=10)

        # Create filter tab
        self.filter_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.filter_frame, text="Filter FASTA")

        # Create merge tab
        self.merge_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.merge_frame, text="Merge FASTA Files")

        # Filter tab variables
        self.input_var = tk.StringVar()
        self.output_var = tk.StringVar()
        self.patterns_var = tk.StringVar(value="##")
        self.regex_var = tk.BooleanVar(value=False)
        self.case_var = tk.BooleanVar(value=False)
        self.report_var = tk.BooleanVar(value=True)

        # Merge tab variables
        self.merge_files = []  # List of Path objects
        self.merge_output_var = tk.StringVar()
        self.dedupe_var = tk.StringVar(value="none")
        self.prefix_var = tk.BooleanVar(value=False)
        self.merge_report_var = tk.BooleanVar(value=True)

        self._build_filter_ui()
        self._build_merge_ui()

    def _build_filter_ui(self):
        pad = {"padx": 10, "pady": 8}

        frm = ttk.Frame(self.filter_frame)
        frm.pack(fill="both", expand=True, **pad)

        # Row 1: Input
        row1 = ttk.Frame(frm)
        row1.pack(fill="x", **pad)
        ttk.Label(row1, text="Input FASTA:").pack(side="left")
        ttk.Entry(row1, textvariable=self.input_var).pack(side="left", fill="x", expand=True, padx=6)
        ttk.Button(row1, text="Browse…", command=self.choose_input).pack(side="left")

        # Row 2: Output
        row2 = ttk.Frame(frm)
        row2.pack(fill="x", **pad)
        ttk.Label(row2, text="Output FASTA:").pack(side="left")
        ttk.Entry(row2, textvariable=self.output_var).pack(side="left", fill="x", expand=True, padx=6)
        ttk.Button(row2, text="Browse…", command=self.choose_output).pack(side="left")

        # Row 3: Patterns
        row3 = ttk.Frame(frm)
        row3.pack(fill="x", **pad)
        ttk.Label(row3, text="Patterns (comma-separated):").pack(side="left")
        ttk.Entry(row3, textvariable=self.patterns_var).pack(side="left", fill="x", expand=True, padx=6)

        # Row 4: Options
        row4 = ttk.Frame(frm)
        row4.pack(fill="x", **pad)
        ttk.Checkbutton(row4, text="Use regular expressions", variable=self.regex_var).pack(side="left")
        ttk.Checkbutton(row4, text="Case sensitive", variable=self.case_var).pack(side="left")
        ttk.Checkbutton(row4, text="Save removal report (.txt)", variable=self.report_var).pack(side="left")

        # Row 5: Actions
        row5 = ttk.Frame(frm)
        row5.pack(fill="x", **pad)
        ttk.Button(row5, text="Run Filter", command=self.run_filter).pack(side="left")

        # Help text
        help_txt = (
            "Tips:\n"
            "• Enter one or more patterns separated by commas. Example: ##,gnl|ECOLI|ABC\n"
            "• If 'Use regular expressions' is ON, each pattern is a regex.\n"
            "• Only headers starting with '>' are checked; sequences are preserved as-is for kept entries.\n"
        )
        row6 = ttk.Frame(frm)
        row6.pack(fill="both", expand=True, **pad)
        ttk.Label(row6, text=help_txt, justify="left").pack(anchor="w")

    def _build_merge_ui(self):
        pad = {"padx": 10, "pady": 8}

        frm = ttk.Frame(self.merge_frame)
        frm.pack(fill="both", expand=True, **pad)

        # Row 1: File list with scrollbar
        row1 = ttk.Frame(frm)
        row1.pack(fill="both", expand=True, **pad)
        ttk.Label(row1, text="Input FASTA files to merge:").pack(anchor="w")
        
        list_frame = ttk.Frame(row1)
        list_frame.pack(fill="both", expand=True, pady=4)
        
        scrollbar = ttk.Scrollbar(list_frame)
        scrollbar.pack(side="right", fill="y")
        
        self.merge_listbox = tk.Listbox(list_frame, yscrollcommand=scrollbar.set, height=6)
        self.merge_listbox.pack(side="left", fill="both", expand=True)
        scrollbar.config(command=self.merge_listbox.yview)
        
        # Buttons for file list management
        btn_frame = ttk.Frame(row1)
        btn_frame.pack(fill="x", pady=4)
        ttk.Button(btn_frame, text="Add Files…", command=self.add_merge_files).pack(side="left", padx=2)
        ttk.Button(btn_frame, text="Remove Selected", command=self.remove_merge_file).pack(side="left", padx=2)
        ttk.Button(btn_frame, text="Clear All", command=self.clear_merge_files).pack(side="left", padx=2)

        # Row 2: Output
        row2 = ttk.Frame(frm)
        row2.pack(fill="x", **pad)
        ttk.Label(row2, text="Output merged FASTA:").pack(side="left")
        ttk.Entry(row2, textvariable=self.merge_output_var).pack(side="left", fill="x", expand=True, padx=6)
        ttk.Button(row2, text="Browse…", command=self.choose_merge_output).pack(side="left")

        # Row 3: Options
        row3 = ttk.Frame(frm)
        row3.pack(fill="x", **pad)
        
        ttk.Label(row3, text="Deduplication:").pack(side="left", padx=(0, 6))
        ttk.Radiobutton(row3, text="None", variable=self.dedupe_var, value="none").pack(side="left")
        ttk.Radiobutton(row3, text="By Header", variable=self.dedupe_var, value="header").pack(side="left")
        ttk.Radiobutton(row3, text="By Sequence", variable=self.dedupe_var, value="sequence").pack(side="left")

        # Row 4: More options
        row4 = ttk.Frame(frm)
        row4.pack(fill="x", **pad)
        ttk.Checkbutton(row4, text="Add file prefix to headers", variable=self.prefix_var).pack(side="left")
        ttk.Checkbutton(row4, text="Save merge report (.txt)", variable=self.merge_report_var).pack(side="left")

        # Row 5: Actions
        row5 = ttk.Frame(frm)
        row5.pack(fill="x", **pad)
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
        row6.pack(fill="x", **pad)
        ttk.Label(row6, text=help_txt, justify="left").pack(anchor="w")

    def choose_input(self):
        path = filedialog.askopenfilename(
            title="Choose FASTA file",
            filetypes=[("FASTA files", "*.fasta *.fa *.faa *.fna *.fas *.fsa"), ("All files", "*.*")],
        )
        if path:
            self.input_var.set(path)
            # Suggest an output next to it
            in_path = Path(path)
            suggested = in_path.with_suffix(in_path.suffix + ".filtered.fasta")
            if not self.output_var.get():
                self.output_var.set(str(suggested))

    def choose_output(self):
        path = filedialog.asksaveasfilename(
            title="Save filtered FASTA as…",
            defaultextension=".fasta",
            filetypes=[("FASTA files", "*.fasta *.fa *.faa *.fna *.fas *.fsa"), ("All files", "*.*")],
        )
        if path:
            self.output_var.set(path)

    def add_merge_files(self):
        paths = filedialog.askopenfilenames(
            title="Choose FASTA files to merge",
            filetypes=[("FASTA files", "*.fasta *.fa *.faa *.fna *.fas *.fsa"), ("All files", "*.*")],
        )
        if paths:
            for path in paths:
                p = Path(path)
                if p not in self.merge_files:
                    self.merge_files.append(p)
                    self.merge_listbox.insert(tk.END, p.name)
            
            # Suggest output if not set
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
            filetypes=[("FASTA files", "*.fasta *.fa *.faa *.fna *.fas *.fsa"), ("All files", "*.*")],
        )
        if path:
            self.merge_output_var.set(path)

    def run_filter(self):
        try:
            input_path = Path(self.input_var.get()).expanduser()
            output_path = Path(self.output_var.get()).expanduser()
            patterns_raw = self.patterns_var.get().strip()
            if not input_path or not input_path.exists():
                messagebox.showerror("Error", "Please choose a valid input FASTA file.")
                return
            if not output_path:
                messagebox.showerror("Error", "Please choose an output FASTA file path.")
                return
            if not patterns_raw:
                messagebox.showerror("Error", "Please enter at least one pattern.")
                return

            patterns = [p.strip() for p in patterns_raw.split(",") if p.strip()]
            report_path = None
            if self.report_var.get():
                report_path = output_path.with_suffix(output_path.suffix + ".removed.txt")

            kept, removed = filter_fasta(
                input_path=input_path,
                output_path=output_path,
                patterns=patterns,
                use_regex=self.regex_var.get(),
                case_sensitive=self.case_var.get(),
                save_report_path=report_path,
            )
            msg = f"Done!\nKept entries: {kept}\nRemoved entries: {removed}"
            if report_path:
                msg += f"\n\nReport saved to:\n{report_path}"
            messagebox.showinfo("FASTA Header Filter", msg)
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def run_merge(self):
        try:
            if not self.merge_files:
                messagebox.showerror("Error", "Please add at least one FASTA file to merge.")
                return
            
            output_path = Path(self.merge_output_var.get()).expanduser()
            if not output_path:
                messagebox.showerror("Error", "Please choose an output FASTA file path.")
                return

            report_path = None
            if self.merge_report_var.get():
                report_path = output_path.with_suffix(output_path.suffix + ".merge_report.txt")

            written, skipped = merge_fasta_files(
                input_paths=self.merge_files,
                output_path=output_path,
                deduplicate=self.dedupe_var.get(),
                add_prefix=self.prefix_var.get(),
                save_report_path=report_path,
            )
            
            msg = f"Done!\nTotal entries written: {written}"
            if skipped > 0:
                msg += f"\nDuplicate entries skipped: {skipped}"
            if report_path:
                msg += f"\n\nReport saved to:\n{report_path}"
            messagebox.showinfo("FASTA Merge", msg)
        except Exception as e:
            messagebox.showerror("Error", str(e))


def main():
    # If run with args, allow headless usage:
    #   python filter_fasta_gui.py input.fasta output.fasta "pat1,pat2" --regex --case
    # Otherwise show the GUI.
    import argparse

    parser = argparse.ArgumentParser(description="Filter FASTA entries by header patterns.")
    parser.add_argument("input", nargs="?", help="Input FASTA")
    parser.add_argument("output", nargs="?", help="Output FASTA")
    parser.add_argument("patterns", nargs="?", help="Comma-separated patterns")
    parser.add_argument("--regex", action="store_true", help="Treat patterns as regular expressions")
    parser.add_argument("--case", action="store_true", help="Case-sensitive matching")
    parser.add_argument("--report", action="store_true", help="Save a removal report next to output")

    args = parser.parse_args()

    if args.input and args.output and args.patterns:
        kept, removed = filter_fasta(
            input_path=Path(args.input),
            output_path=Path(args.output),
            patterns=[p.strip() for p in args.patterns.split(",") if p.strip()],
            use_regex=args.regex,
            case_sensitive=args.case,
            save_report_path=Path(args.output).with_suffix(Path(args.output).suffix + ".removed.txt") if args.report else None,
        )
        print(f"Kept entries: {kept}")
        print(f"Removed entries: {removed}")
    else:
        app = App()
        app.mainloop()


if __name__ == "__main__":
    main()
