# Fixed inputs

The six files in `fixed/` are the recorded common RE/DE inputs from the model
software accompanying Bordalo, Gennaioli, Shleifer, and Terry (2026),
*Real Credit Cycles*, https://doi.org/10.3886/E237827V1.

- `Asimshocks.txt`: aggregate simulation draws.
- `zfirmshocks.txt`: firm simulation draws.
- `upifirmshocks1.txt` and `upifirmshocks2.txt`: measurement-noise simulation inputs.
- `_GMMdatamom.txt`: target moments.
- `_GMMdatawgt.txt`: moment weight matrix.

These are simulation and objective inputs, not original analyst-level records.
Their SHA-256 identities are in `results/integrated-manifest.json`. They are
identical in the recorded RE and DE source folders. See `NOTICE` for attribution.
