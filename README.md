# Documentation workflow

The tutorial source is `example.qmd`. Quarto is used as the generator, and the rendered HTML is committed so CI can deploy documentation without executing Julia notebooks.

To regenerate the committed documentation after editing `example.qmd`, run from the repository root:

```sh
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
quarto render docs/example.qmd
```

This updates `docs/index.html` and any files under `docs/example_files/`. Commit those generated files together with the `.qmd` source.

Longer-form package documentation can later be assembled with Documenter.jl, while keeping Quarto as the tutorial generator.
