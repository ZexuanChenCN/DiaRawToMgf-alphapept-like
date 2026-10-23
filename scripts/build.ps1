[CmdletBinding()]
param([switch]$CetCompatibility)
Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'
$repoRoot = Split-Path -Parent $PSScriptRoot
$previousCompilerHost = $env:DIA_CSC_DOTNET
$previousCompilerDll = $env:DIA_CSC_DLL
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
    $publishArgs = @('publish', 'DiaRawToMgf.csproj', '-c', 'Release', '--no-restore', '--no-self-contained', '-o', 'artifacts/publish')
    if ($CetCompatibility) {
        $sdkLine = & dotnet --list-sdks | Where-Object { $_ -match ('^' + [regex]::Escape($sdkVersion) + '\s+\[') } | Select-Object -First 1
        if (-not $sdkLine) { throw "Cannot locate installed SDK $sdkVersion." }
        $sdkMatch = [regex]::Match($sdkLine, '\[(.+)\]\s*$')
        $env:DIA_CSC_DOTNET = (Get-Command dotnet -CommandType Application | Select-Object -First 1).Source
        $env:DIA_CSC_DLL = Join-Path $sdkMatch.Groups[1].Value "$sdkVersion\Roslyn\bincore\csc.dll"
        if (-not (Test-Path -LiteralPath $env:DIA_CSC_DLL)) { throw "Compiler DLL not found: $env:DIA_CSC_DLL" }
        Write-Host 'CET compatibility: use dotnet csc.dll; disable CET marking only for this application.'
        $publishArgs += @('-p:UseSharedCompilation=false', "-p:CscToolPath=$PSScriptRoot", '-p:CscToolExe=csc-via-dotnet.cmd', '-p:CETCompat=false')
    }
    & dotnet @publishArgs
    if ($LASTEXITCODE -ne 0) { throw 'Publish failed. For a CET startup error, rerun this script with -CetCompatibility; see README.' }
    & dotnet artifacts/publish/DiaRawToMgf.dll --self-test
    if ($LASTEXITCODE -ne 0) { throw 'Self-test failed.' }
    Write-Host "Ready: $(Join-Path $repoRoot 'artifacts/publish/DiaRawToMgf.exe')"
}
finally {
    $env:DIA_CSC_DOTNET = $previousCompilerHost
    $env:DIA_CSC_DLL = $previousCompilerDll
    Pop-Location
}
