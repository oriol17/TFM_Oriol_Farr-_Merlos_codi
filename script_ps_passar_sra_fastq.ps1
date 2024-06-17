$output_dir = "D:\Master\TFM\SRAs"
$srr_list = @("SRR18283143", "SRR18283144", "SRR18283145", "SRR18283146", "SRR18283147", "SRR18283148", "SRR18283149", "SRR18283150", "SRR18283151", "SRR18283152", "SRR18283153", "SRR18283154", "SRR18283155", "SRR18283156", "SRR18283157", "SRR18283158", "SRR18283159", "SRR18283160", "SRR18283161", "SRR18283162", "SRR18283163", "SRR18283164", "SRR18283165", "SRR18283166", "SRR18283167", "SRR18283168", "SRR18283169", "SRR18283170", "SRR18283171")

cd C:\Users\oriol\Downloads\sratoolkit.3.1.1-win64\sratoolkit.3.1.1-win64\bin

foreach ($srr in $srr_list) {
    Write-Output "Processing $srr"
    .\fasterq-dump $srr -O $output_dir
}

Write-Output "All downloads completed."