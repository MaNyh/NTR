===
CascadeCase_gwk_rans_trace_parastud
===

This example can be used to reproduce a RANS parameterstudy in TRACE

Usage
-------------

- define parameters in settings.yml
- run run_create_geometry
- run run_meshing
- run run_parastudcreator
- copy the cgns mesh into the input-dir of the simulation
- run gmcplay with the journal
- run the simulation

Features
-------------

The shown methods be used for parameterstudies that are defined with ascii-files.
This makes it suitable for any kind of numerical solver
