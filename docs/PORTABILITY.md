# Portability notes

The repository contains 54 R analysis scripts. The [manifest](../script_manifest.csv) records their repository paths and SHA-256 checksums.

[Path configuration](../config/paths.example.R) defaults to the repository's `data/` and `results/` directories. Set `OAGLY_DATA_ROOT` or `OAGLY_RESULTS_ROOT` to use another location. You can also copy the example to `config/paths.R` for local configuration. Git ignores this file. When launching a script outside the repository root, set `OAGLY_PROJECT_ROOT` to the repository directory. See the [workflow](WORKFLOW.md#rds-handoffs) for input locations and intermediate objects.

The [dependency list](../requirements.R) is a starting point. Package versions vary across historical analyses. The complete workflow has yet to be rerun end to end in a clean environment with the original inputs. The local verification tools check R syntax and scan the repository for bulky data files, credentials and local path references.

```powershell
./tools/check_r_syntax.ps1
./tools/check_repository.ps1
```
