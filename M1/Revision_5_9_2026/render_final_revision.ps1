$ErrorActionPreference = 'Stop'

$root = Split-Path -Parent $MyInvocation.MyCommand.Path
$jobs = @(
    @{ Source = 'Manuscript- Water and Ecology_Fully_Revised_Clean.docx'; Output = '_qa_final3_main_clean' },
    @{ Source = 'Manuscript- Water and Ecology_Fully_Revised_Colour_Marked.docx'; Output = '_qa_final3_main_marked' },
    @{ Source = 'SupplementaryInformation_Fully_Revised_Clean.docx'; Output = '_qa_final3_supp_clean' },
    @{ Source = 'SupplementaryInformation_Fully_Revised_Colour_Marked.docx'; Output = '_qa_final3_supp_marked' },
    @{ Source = 'Response_to_Reviewers_Fully_Addressed.docx'; Output = '_qa_final3_response' }
)

$word = New-Object -ComObject Word.Application
$word.Visible = $false
$word.DisplayAlerts = 0
try {
    foreach ($job in $jobs) {
        $source = Join-Path $root $job.Source
        $outputDirectory = Join-Path $root $job.Output
        New-Item -ItemType Directory -Path $outputDirectory -Force | Out-Null
        Get-ChildItem -LiteralPath $outputDirectory -File | Where-Object {
            $_.Name -like 'page-*.png' -or
            $_.Name -like 'contact-*.png' -or
            $_.Name -eq 'render.pdf'
        } | Remove-Item -Force
        $pdf = Join-Path $outputDirectory 'render.pdf'
        $document = $word.Documents.Open($source, $false, $true)
        try {
            $document.ExportAsFixedFormat($pdf, 17)
        }
        finally {
            $document.Close($false)
        }
        & pdftoppm -png -r 150 $pdf (Join-Path $outputDirectory 'page')
        if ($LASTEXITCODE -ne 0) {
            throw "pdftoppm failed for $source"
        }
    }
}
finally {
    $word.Quit()
    [System.Runtime.InteropServices.Marshal]::FinalReleaseComObject($word) | Out-Null
}

$jobs | ForEach-Object {
    $directory = Join-Path $root $_.Output
    [pscustomobject]@{
        Set = $_.Output
        Pages = @(Get-ChildItem -LiteralPath $directory -Filter 'page-*.png').Count
        PdfBytes = (Get-Item -LiteralPath (Join-Path $directory 'render.pdf')).Length
    }
}
