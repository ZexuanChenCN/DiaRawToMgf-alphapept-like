# Downloads the exact official SDK packages used for the validated RAW run.
[CmdletBinding()]
param()
Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'
[Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12
$repoRoot = Split-Path -Parent $PSScriptRoot
$feed = Join-Path $repoRoot '.nuget/vendor'
$revision = '345aafcb46c944b843ca49880b5644bff6285f61'
$baseUri = "https://raw.githubusercontent.com/thermofisherlsms/RawFileReader/$revision/Libs/NetCore/Net8"
$packages = @(
    @{ Name = 'ThermoFisher.CommonCore.Data.8.0.6.nupkg'; Sha256 = 'AA6C6A170A43200F8C2DE3756C94DF5D20CCC8212B4D64F35018D3E2B69B6521' },
    @{ Name = 'ThermoFisher.CommonCore.RawfileReader.8.0.6.nupkg'; Sha256 = '61921FACFC165239B4F7752BD5DF01E8FE8764067F24DEFE1F06A094FCD208DC' }
)
New-Item -ItemType Directory -Force -Path $feed | Out-Null
foreach ($package in $packages) {
    $path = Join-Path $feed $package.Name
    if (Test-Path -LiteralPath $path) {
        if ((Get-FileHash -LiteralPath $path -Algorithm SHA256).Hash -ne $package.Sha256) {
            throw "Existing SDK package has an unexpected SHA256: $path"
        }
        Write-Host "Verified: $($package.Name)"
        continue
    }
    $partial = "$path.partial"
    Write-Host "Downloading: $($package.Name)"
    Invoke-WebRequest -UseBasicParsing -Uri "$baseUri/$($package.Name)" -OutFile $partial
    if ((Get-FileHash -LiteralPath $partial -Algorithm SHA256).Hash -ne $package.Sha256) {
        throw "SDK download failed SHA256 verification: $partial"
    }
    Move-Item -LiteralPath $partial -Destination $path
}
Write-Host 'Thermo SDK 8.0.6 ready. Vendor license: https://github.com/thermofisherlsms/RawFileReader/blob/345aafcb46c944b843ca49880b5644bff6285f61/License.doc'
