param([string]$Python = 'python', [string]$ExtraPackages = '')
$ErrorActionPreference = 'Stop'
Push-Location -LiteralPath $PSScriptRoot
try {
    foreach ($figureScript in @('make_charts_v3.py','make_singlecol_figures.py','make_geometry_v3.py','make_limitations_v3.py','make_advective_figure_v3.py','compose_figures_v3.py')) {
        $runnerArguments = @('analysis/run_python.py')
        if ($ExtraPackages) { $runnerArguments += @('--extra-packages',$ExtraPackages) }
        $runnerArguments += "analysis/$figureScript"
        if ($figureScript -eq 'make_singlecol_figures.py') { $runnerArguments += '.' }
        & $Python @runnerArguments
        if ($LASTEXITCODE -ne 0) { throw "Figure rebuild failed: $figureScript" }
    }
} finally { Pop-Location }
