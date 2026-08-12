"""
===========================================================================
 TractQuant — "Validate" page integration patch
===========================================================================
Put `section_matcher.py` next to TractQuant.pyw and pipeline.py, then make the
four small edits below inside the TractQuant class in TractQuant.pyw.

Edit 1, 2, 3 register the page. Edit 4 makes the actual analysis read exactly
the image pairing you confirmed on the Validate page (so the run can never
silently disagree with what you eyeballed).
===========================================================================
"""

# ───────────────────────────────────────────────────────────────────────
# EDIT 1 — add the page to the PAGES list (put it right after "data")
#
#   PAGES=[("setup","⚙  Setup"),("data","📂  Load data"),
#          ("validate","✅  Validate"),          # <-- ADD THIS LINE
#          ("config","🔧  Configure"),("run","▶  Run"),
#          ("results","📊  Results"),("export","💾  Export")]
# ───────────────────────────────────────────────────────────────────────


# ───────────────────────────────────────────────────────────────────────
# EDIT 2 — add a title for it in _show()'s `titles` dict:
#
#   "validate":("Validate alignment","Confirm each section is paired with the right image"),
# ───────────────────────────────────────────────────────────────────────


# ───────────────────────────────────────────────────────────────────────
# EDIT 3 — add it to the dispatch dict at the bottom of _show():
#
#   {"setup":self._p_setup,"data":self._p_data,"validate":self._p_validate,
#    "config":self._p_config,"run":self._p_run,
#    "results":self._p_results,"export":self._p_export}[pid]()
# ───────────────────────────────────────────────────────────────────────


# ───────────────────────────────────────────────────────────────────────
# EDIT 4 — paste these two methods into the TractQuant class
#          (anywhere among the other _p_* methods)
# ───────────────────────────────────────────────────────────────────────

def _validate_all(self):
    """Run strict matching for every animal. Returns {animal_id: records}."""
    import importlib, pipeline as pl, section_matcher as sm
    importlib.reload(pl); importlib.reload(sm)
    out = {}
    for an in self.animals:
        sections = pl.parse_visual_align(Path(an["json_path"]))
        recs = sm.match_sections_to_images(
            sections, an["ch1_dir"], an["ch2_dir"], animal_id=an["id"]
        )
        out[an["id"]] = recs
    self._match_cache = out          # reused by the run loop (Edit 4b)
    return out

def _p_validate(self):
    import section_matcher as sm
    p = self.pf
    if not self.animals:
        lbl(p, "Add animals first (Load data).", color=FG2, bg=BG2).pack(pady=40)
        return

    intro = card(p); intro.pack(fill="x", padx=16, pady=(14, 6))
    sec_title(intro, "Section ↔ image check")
    lbl(intro,
        "Each VisualAlign section is matched to its Ch1/Ch2 image by section id "
        "(animal prefix stripped, so digits in the animal name can't collide). "
        "Build the overlay sheet to see the CCF atlas warped onto each real "
        "section — confirm the outlines hug the anatomy before running.",
        color=FG2, bg=BG2, wraplength=720).pack(anchor="w", padx=12, pady=(0, 10))

    status_card = card(p); status_card.pack(fill="x", padx=16, pady=6)
    sec_title(status_card, "Per-animal result")
    rows = tk.Frame(status_card, bg=BG2); rows.pack(fill="x", padx=12, pady=(0, 10))

    note = tk.Label(p, text="", font=("Segoe UI", 10), bg=BG2, fg=FG2)
    note.pack(anchor="w", padx=16)

    def render(results):
        for w in rows.winfo_children(): w.destroy()
        any_problem = False
        for an in self.animals:
            recs = results.get(an["id"], [])
            s = sm.summarize(recs)
            problem = s["problems"] > 0
            any_problem = any_problem or problem
            r = tk.Frame(rows, bg=BG3, highlightthickness=1,
                         highlightbackground=(RED if problem else BORDER))
            r.pack(fill="x", pady=3)
            lbl(r, an["id"], bold=True, bg=BG3).pack(side="left", padx=10, pady=8)
            txt = (f"{s['ok']}/{s['total']} matched"
                   if not problem else
                   f"{s['ok']}/{s['total']} matched · "
                   + " · ".join(f"{s[k]} {k}" for k in
                                ("missing", "duplicate", "ambiguous", "channel_mismatch")
                                if s.get(k)))
            tk.Label(r, text=txt, font=("Consolas", 9), bg=BG3,
                     fg=(GREEN if not problem else AMBER)).pack(side="left", padx=6)
        note.config(
            text=("✓ All animals matched cleanly — safe to run."
                  if not any_problem else
                  "⚠ Some sections need attention. Open a contact sheet to inspect."),
            fg=(GREEN if not any_problem else AMBER))

    def build_sheets():
        gen_btn.config(state="disabled", text="Loading atlas…")
        def _go():
            import importlib, webbrowser
            import pipeline as pl, atlas_manager as am, section_matcher as sm
            importlib.reload(pl); importlib.reload(sm)
            try:
                # atlas loaded ONCE (same as the run does)
                atlas = pl.CCFAtlas(str(am.get_atlas_volume()), str(am.get_region_csv()))
                resize_map = {"25%": 0.25, "50%": 0.5, "100%": 1.0}
                resize = resize_map.get(self.v_resize.get().split("—")[0].strip(), 0.25)
                outdir = (APP_DIR / "results" /
                          (self.v_project.get().replace(" ", "_") or "project"))
                qc_dir = outdir / "alignment_qc"
                results, first = {}, None
                for an in self.animals:
                    sections = pl.parse_visual_align(Path(an["json_path"]))
                    recs = sm.match_sections_to_images(
                        sections, an["ch1_dir"], an["ch2_dir"], animal_id=an["id"])
                    def prog(i, n, aid=an["id"]):
                        self.after(0, lambda: gen_btn.config(
                            text=f"{aid}: overlay {i}/{n}…"))
                    recs = sm.build_overlay_records(
                        recs, sections, atlas, outdir, an["id"], pl,
                        resize_factor=resize, thumb_px=440, progress=prog)
                    results[an["id"]] = recs
                    path = qc_dir / f"{an['id']}_validation.html"
                    sm.build_contact_sheet(recs, path, animal_id=an["id"], thumb_px=440)
                    first = first or path
                self._match_cache = results
                self.after(0, lambda: render(results))
                if first:
                    webbrowser.open(first.resolve().as_uri())
            except Exception as e:
                import traceback; tb = traceback.format_exc()
                self.after(0, lambda: note.config(text=f"Error: {e}", fg=RED))
                print(tb)
            finally:
                self.after(0, lambda: gen_btn.config(
                    state="normal", text="Build & open overlay sheet  ↗"))
        threading.Thread(target=_go, daemon=True).start()

    bf = tk.Frame(p, bg=BG2); bf.pack(fill="x", padx=16, pady=12)
    action_btn(bf, "Re-check matches", lambda: render(self._validate_all()),
               color=PURP2).pack(side="left", padx=(0, 8))
    gen_btn = action_btn(bf, "Build & open overlay sheet  ↗", build_sheets, color=GREEN)
    gen_btn.pack(side="left")
    lbl(bf, "  (overlays the VisualAlign atlas on each section to confirm the alignment)",
        color=FG2, bg=BG2).pack(side="left", padx=6)

    render(self._validate_all())     # initial fast filename check (no atlas/images)


# ───────────────────────────────────────────────────────────────────────
# EDIT 4b — make the RUN use the validated pairing.
#
# In _run_thread(), replace the section loop. Current code:
#
#     sections=pl.parse_visual_align(Path(an["json_path"]))
#     section_rows=[]
#     for filename,sec in sections.items():
#         nr=sec["nr"]
#         try:
#             ch1,ch2=pl.load_section_images(Path(an["ch1_dir"]),Path(an["ch2_dir"]),nr,filename)
#             ...
#
# Replace with:
#
#     import section_matcher as sm
#     sections=pl.parse_visual_align(Path(an["json_path"]))
#     recs={r["json_file"]: r for r in
#           sm.match_sections_to_images(sections, an["ch1_dir"], an["ch2_dir"],
#                                       animal_id=an["id"])}
#     section_rows=[]
#     for filename,sec in sections.items():
#         nr=sec["nr"]; m=recs[filename]
#         if m["status"] in ("missing","duplicate","channel_mismatch"):
#             self._log(f"  s{nr:03d} ⚠ skipped — {m['status']} ({m['notes']})",AMBER)
#             continue
#         try:
#             ch1,ch2=sm.load_validated_pair(m)        # reads the confirmed files
#             sdf=pl.measure_section(sec,ch1,ch2,atlas,resize)
#             qc=pl.save_alignment_overlay(sec,ch1,ch2,atlas,output_dir,an["id"],nr,filename,resize)
#             sdf["section_nr"]=nr; sdf["section_file"]=filename
#             section_rows.append(sdf)
#             ... (keep the existing logging lines)
#         except Exception as e:
#             self._log(f"  s{nr:03d} ✗ {e}",RED)
#
# This way the numbers are computed from exactly the images you validated.
# ───────────────────────────────────────────────────────────────────────
