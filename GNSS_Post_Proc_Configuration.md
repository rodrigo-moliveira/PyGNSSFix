# PositionFix Configuration Reference

This document describes the fields available in the PositionFix JSON configuration files (`config.json`, `config_post.json`).  
Each section lists the parameters, their purpose, and possible values.

---

## General

- `mode`: Positioning mode to use.  
  Possible values: `"SPP"`, `"PR-PPP"`, `"CP-PPP"`.

- `output_dir`: Path to the directory where results will be saved.

---

## Input Files

- `obs_file`: Path to the RINEX Observation file (e.g., `"station.obs"`).  
- `nav_file`: Path to the RINEX Navigation file (e.g., `"brdc.nav"`).  
- `clk_file`: Path to the RINEX CLK file containing precise satellite clocks.  
- `sp3_file`: Path to the SP3 file containing precise satellite orbits.  
- `ionex_file`: Path to the IONEX file with global ionosphere maps.  
- `antex_file`: Path to ANTEX file for antenna phase center offsets/variations (optional).  

---