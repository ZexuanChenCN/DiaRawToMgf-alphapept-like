[CmdletBinding()]
param()
Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'
$repoRoot = Split-Path -Parent $PSScriptRoot
Push-Location $repoRoot
try {
    $sdkVersion = & dotnet --version
    if ($LASTEXITCODE -ne 0 -or $sdkVersion -notmatch '^10\.') {
        throw 'Install the .NET 10 SDK, reopen PowerShell, and try again.'
    }
    & (Join-Path $PSScriptRoot 'prepare-sdk.ps1')
    $restoreArgs = @('restore', 'DiaRawToMgf.csproj', '--configfile', 'NuGet.Config')
    if (Test-Path -LiteralPath 'packages.lock.json') { $restoreArgs += '--locked-mode' }
    & dotnet @restoreArgs
    if ($LASTEXITCODE -ne 0) { throw 'NuGet restore failed.' }
    & dotnet publish DiaRawToMgf.csproj -c Release --no-restore --no-self-contained -o artifacts/publish
    if ($LASTEXITCODE -ne 0) { throw 'Publish failed.' }
    & dotnet artifacts/publish/DiaRawToMgf.dll --self-test
    if ($LASTEXITCODE -ne 0) { throw 'Self-test failed.' }
    Write-Host "Ready: $(Join-Path $repoRoot 'artifacts/publish/DiaRawToMgf.exe')"
}
finally { Pop-Location }
