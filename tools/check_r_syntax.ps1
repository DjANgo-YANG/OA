param(
    [string]$RepositoryRoot = (Split-Path -Parent $PSScriptRoot)
)

$ErrorActionPreference = 'Stop'
$root = [System.IO.Path]::GetFullPath($RepositoryRoot)
$scriptsRoot = Join-Path $root 'scripts'
$rscript = (Get-Command Rscript -ErrorAction Stop).Source
$files = @(Get-ChildItem -LiteralPath $scriptsRoot -Recurse -File -Filter '*.R')
$failures = @()

foreach ($file in $files) {
    $content = Get-Content -LiteralPath $file.FullName -Raw -Encoding UTF8
    $startInfo = New-Object System.Diagnostics.ProcessStartInfo
    $startInfo.FileName = $rscript
    $startInfo.Arguments = '--vanilla -e "parse(file(''stdin''), keep.source=TRUE)"'
    $startInfo.UseShellExecute = $false
    $startInfo.CreateNoWindow = $true
    $startInfo.RedirectStandardInput = $true
    $startInfo.RedirectStandardOutput = $true
    $startInfo.RedirectStandardError = $true
    $process = New-Object System.Diagnostics.Process
    $process.StartInfo = $startInfo
    [void]$process.Start()
    $process.StandardInput.Write($content)
    $process.StandardInput.Close()
    [void]$process.StandardOutput.ReadToEnd()
    [void]$process.StandardError.ReadToEnd()
    $process.WaitForExit()
    if ($process.ExitCode -ne 0) {
        $failures += $file.FullName.Substring($root.Length + 1)
    }
    $process.Dispose()
}

Write-Host "Parsed $($files.Count) R scripts."
if ($failures.Count) {
    Write-Host "Syntax failures: $($failures.Count)"
    $failures | ForEach-Object { Write-Host "- $_" }
    exit 1
}

Write-Host 'All R scripts passed parse-only syntax validation.'
exit 0
