param(
    [string]$RepositoryRoot = (Split-Path -Parent $PSScriptRoot),
    [int]$MaxFileSizeMB = 50
)

$ErrorActionPreference = 'Stop'
$root = [System.IO.Path]::GetFullPath($RepositoryRoot)
$maxBytes = $MaxFileSizeMB * 1MB

Write-Host "Repository: $root"

$gitFiles = @(& git -C $root ls-files --cached --others --exclude-standard 2>$null)
if ($LASTEXITCODE -eq 0 -and $gitFiles.Count) {
    $allFiles = @($gitFiles | ForEach-Object {
        Get-Item -LiteralPath (Join-Path $root $_)
    })
} else {
    $allFiles = @(Get-ChildItem -LiteralPath $root -Recurse -File -Force |
        Where-Object { $_.FullName -notmatch '[\\/]\.git[\\/]' })
}

$large = @($allFiles | Where-Object { $_.Length -gt $maxBytes })
$blockedExtensions = @(
    '.rds', '.rda', '.rdata', '.loom', '.h5', '.h5ad', '.hdf5',
    '.bam', '.bai', '.sam', '.cram', '.fastq', '.fq', '.sra',
    '.ai', '.psd', '.doc', '.docx'
)
$scientificData = @($allFiles | Where-Object {
    $name = $_.Name.ToLowerInvariant()
    $_.Extension.ToLowerInvariant() -in $blockedExtensions -or
    $name.EndsWith('.fastq.gz') -or $name.EndsWith('.fq.gz')
})

$textExtensions = @('.r', '.rmd', '.qmd', '.sh', '.md', '.txt', '.yml', '.yaml', '.json', '.csv', '.ps1')
$textFiles = @($allFiles | Where-Object { $_.Extension.ToLowerInvariant() -in $textExtensions })

$secretPattern = '(?i)(api[_-]?key|access[_-]?token|client[_-]?secret|password|passwd|bearer)\s*[:=]\s*["''][^"'']{6,}'
$pathPattern = '(?i)(?:^|["''\s=])(?:[A-Z]:[\\/]|/home/|/root/|/mnt/|/media/)'

$secretHits = @()
$pathHits = @()
foreach ($file in $textFiles) {
    $lineNumber = 0
    foreach ($line in [System.IO.File]::ReadLines($file.FullName)) {
        $lineNumber++
        if ($line -match $secretPattern) {
            $secretHits += "$($file.FullName.Substring($root.Length + 1)):$lineNumber"
        }
        if ($line -match $pathPattern) {
            $pathHits += "$($file.FullName.Substring($root.Length + 1)):$lineNumber"
        }
    }
}

Write-Host "Files scanned: $($allFiles.Count)"
Write-Host "Files over ${MaxFileSizeMB} MB: $($large.Count)"
Write-Host "Blocked scientific/data files: $($scientificData.Count)"
Write-Host "Potential secret assignments: $($secretHits.Count)"
Write-Host "Environment-specific path lines: $($pathHits.Count)"

if ($large.Count) {
    Write-Host "`nLarge files:"
    $large | ForEach-Object { Write-Host "- $($_.FullName.Substring($root.Length + 1))" }
}
if ($scientificData.Count) {
    Write-Host "`nData files requiring review:"
    $scientificData | ForEach-Object { Write-Host "- $($_.FullName.Substring($root.Length + 1))" }
}
if ($secretHits.Count) {
    Write-Host "`nPotential secrets:"
    $secretHits | ForEach-Object { Write-Host "- $_" }
}
if ($pathHits.Count) {
    Write-Host "`nPortability warnings (first 100):"
    $pathHits | Select-Object -First 100 | ForEach-Object { Write-Host "- $_" }
}

if ($large.Count -or $scientificData.Count -or $secretHits.Count) {
    exit 1
}

exit 0
