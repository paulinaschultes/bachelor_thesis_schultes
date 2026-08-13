inputDir = getDirectory("Choose input folder with PNGs");


lutFile = File.openDialog("Choose LUT file (.lut)");
outputDir = getDirectory("Choose output folder for processed images");
list = getFileList(inputDir);
open(lutFile);
getLut(redLUT, greenLUT, blueLUT);
close(); 
for (i = 0; i < list.length; i++) {
    if (endsWith(list[i], ".png")) {
        open(inputDir + list[i]);
        setLut(redLUT, greenLUT, blueLUT);
        run("RGB Color");
        saveAs("PNG", outputDir + list[i]);
        close();
    }
}

print("Linus saved your time <3");