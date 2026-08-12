"""
Background diagnostic for TractQuant.

Double-click Background-Test.command, or run:
    python diagnose_background.py

You can also use command-line arguments:
    python diagnose_background.py <image> <visualign.json> <section_nr>
"""

from pathlib import Path
import sys

from pipeline import diagnose_background_method


def _pick_files_with_dialog():
    import tkinter as tk
    from tkinter import filedialog, simpledialog, messagebox

    root = tk.Tk()
    root.withdraw()
    root.update()

    image_path = filedialog.askopenfilename(
        title="Choose the Ch1 image for one section",
        filetypes=[
            ("Image files", "*.tif *.tiff *.png *.jpg *.jpeg"),
            ("All files", "*.*"),
        ],
    )
    if not image_path:
        return None

    json_path = filedialog.askopenfilename(
        title="Choose the matching VisualAlign JSON",
        filetypes=[("JSON files", "*.json"), ("All files", "*.*")],
    )
    if not json_path:
        return None

    section_nr = simpledialog.askinteger(
        "Section number",
        "Which section number should be tested?",
        minvalue=0,
    )
    if section_nr is None:
        return None

    root.destroy()
    return Path(image_path), Path(json_path), int(section_nr)


def main():
    if len(sys.argv) >= 4:
        image_path = Path(sys.argv[1])
        json_path = Path(sys.argv[2])
        section_nr = int(sys.argv[3])
    else:
        picked = _pick_files_with_dialog()
        if picked is None:
            print("Cancelled.")
            return
        image_path, json_path, section_nr = picked

    print()
    print("Running TractQuant background diagnostic...")
    print(f"Image: {image_path}")
    print(f"VisualAlign JSON: {json_path}")
    print(f"Section: {section_nr}")
    print()

    diagnose_background_method(image_path, json_path, section_nr)

    print()
    print("Done. You can copy these values into your validation notes.")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print()
        print("Background diagnostic failed:")
        print(exc)
        print()
        print("Tip: make sure the section number exists in the VisualAlign JSON")
        print("and that the selected image belongs to that same section.")
        raise
