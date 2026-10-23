[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)][string]$PythonMgf,
    [Parameter(Mandatory = $true)][string]$CSharpMgf
)
Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'
$reference = Get-Item -LiteralPath $PythonMgf
$candidate = Get-Item -LiteralPath $CSharpMgf
if ($reference.FullName -eq $candidate.FullName) { throw 'Supply two different files.' }
$refHash = (Get-FileHash -LiteralPath $reference.FullName -Algorithm SHA256).Hash
$candidateHash = (Get-FileHash -LiteralPath $candidate.FullName -Algorithm SHA256).Hash
$same = $reference.Length -eq $candidate.Length -and $refHash -eq $candidateHash
[pscustomobject]@{
    PythonFile = $reference.FullName
    CSharpFile = $candidate.FullName
    PythonBytes = $reference.Length
    CSharpBytes = $candidate.Length
    PythonSha256 = $refHash
    CSharpSha256 = $candidateHash
    IdenticalLengthAndSha256 = $same
}
if (-not $same) { throw 'MGF size or SHA256 differs. Use qa/compare_csharp_dia.py for a semantic comparison.' }
