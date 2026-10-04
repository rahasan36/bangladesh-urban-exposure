# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

import os, sys, subprocess, time

PKG = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable

STAGES = [
    ("pipeline", ["00b_prepare_vectors.py", "00_synthetic_tests.py", "01_harmonize.py", "02_areas_transitions.py", "03_decomposition.py",
                  "04_pa_kba_buffers.py", "05a_species_prep.py", "05b_species_zonal.py", "05c_species_models.py", "06_hotspots.py",
                  "07_species_special.py", "08_tile_test.py", "09_validation.py"]),
    ("multires", ["01_ladder_glc.py", "02_ladder_glc_strips.py", "03_ladder_ghsl.py", "04_ghsl_thresholds.py", "05_zones_species.py",
                  "06_ghsl_full.py", "07_sensitivities.py", "08_sundarbans_buffers.py"]),
    ("sdm", ["00_build_sdm_grid.py", "01_species_urban_share.py", "02_species_tests_models.py", "03_hotspot_masks.py",
             "04_hotspot_overlaps.py", "05_validation_hotspots.py"]),
    ("analysis", ["hotspot_domain.py", "group_gain_tests.py", "family_mixed_model.py", "share2020_models.py", "matched_rankings.py",
                  "site_reconcile.py", "data_s1.py"]),
    ("figures", ["fig2_data.py", "fig2_distance_profile.py"]),
    ("tables", ["build_si_tables.py"]),
    ("figures", ["fig1_resolution.py", "fig2_pa_kba.py", "fig3_species_status.py", "fig4_species_groups.py", "fig5_hotspots.py",
                 "figS1_validation.py", "figS2_site_ranking.py", "figS3_sundarbans.py", "figS4_S5_S6.py", "figS7_hotspots_products.py"]),
]

def run(script, cwd, args=()):
    for _ in range(1000):
        t0 = time.time()
        p = subprocess.run([PY, script, *args], cwd=cwd, text=True, capture_output=True)
        tail = (p.stdout + p.stderr).strip().splitlines()
        print(f"  {script}: exit {p.returncode}, {time.time() - t0:.0f}s")
        if p.returncode != 0:
            print("\n".join(tail[-30:])); raise SystemExit(f"{script} failed")
        if "CONTINUE" not in p.stdout: return
    raise SystemExit(f"{script} did not finish")

if __name__ == "__main__":
    a = sys.argv[1:] or ["all"]
    if a[0] == "--check":
        run("00a_download_inputs.py", os.path.join(PKG, "pipeline"), ["--check"]); sys.exit()
    only = a[1] if len(a) > 1 else None
    for folder, scripts in STAGES:
        if a[0] not in ("all", folder): continue
        for s in scripts:
            if only and not s.startswith(only): continue
            print(f"[{folder}] {s}"); run(s, os.path.join(PKG, folder))
    print("done")
