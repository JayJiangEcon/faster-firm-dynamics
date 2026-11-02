# CUDA computation and search

Version 0.2.0 provides a complete batch-oriented search pipeline:

- [Build, evaluation, search and SSH execution](cuda-pipeline.md)
- [Version 0.2.0 verification, timing scope and execution identities](cuda-validation-v020.md)
- [Historical version 0.1.0 GPU measurements](gpu-v010.md)

The new pipeline does not need a CPU solution before evaluating a new
candidate. CPU input construction and the CUDA context remain resident across
batches; each candidate starts a cold model solve. Mac numerical finishing preserves the macro-statistic and direct-
objective arithmetic. Measurement-noise profiling
then supplies the objective for a deferred population search.

## Overall speed relative to the original M3 baseline

The large throughput ratios previously reported compare original M3 evaluation
rates with repeated-vector GPU validation batches using pre-exported inputs.
See the [historical table](gpu-v010.md#overall-speed-relative-to-the-original-m3-baseline)
for those observations. They are not version 0.2.0 end-to-end search speedups.
The current validation report records the actual scope of each new timing.
