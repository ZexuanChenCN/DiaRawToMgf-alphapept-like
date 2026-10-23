Validation reports
==================
alphapept_port_* reports refer to the AlphaPept 0.5.3 numerical source port.
The slice report compares all 14 columns generated directly from RAW, RT 10-12
minutes, against an independent Python/AlphaPept oracle. Stage snapshots also
compare hill arrays, bootstrap values, isotope groups and centroid assignments.
Development scripts used for comparison are in ../qa/.
alphapept_port_full_features.json compares all 223,027 features and 14 columns.
alphapept_port_full_mgf.json compares all candidate pairs, charges, ranks,
TITLEs, complete peak arrays and 5,660,503 MS1 curve points; final PASS.
alphapept_port_full_summary.json records the actual pure C# full RAW run:
247.72 seconds total, including 92.71 seconds feature finding, four threads.

Historical reports (retired approximate engine, NOT the current port):
native_comparison.json and native_summary.json use native-hill-isotope-v1.
Their 44.53% candidate agreement motivated replacing that implementation.
association_comparison.json bypasses feature finding using a reference TSV
cache; it validates only association/export, not the source port.
