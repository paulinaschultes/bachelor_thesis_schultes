input = getDirectory("Choose Input Folder");
output = getDirectory("Choose Output Folder");
list = getFileList(input);

for (i = 0; i < list.length; i++) {
    if (endsWith(list[i], ".tif") || endsWith(list[i], ".jpg") || endsWith(list[i], ".png")) {

        open(input + list[i]);

       

        // Set brightness/contrast manually: min=0, max=300
        setMinAndMax(512, 21248);

        // Bake the brightness/contrast into pixel values
        run("Apply LUT");

        saveAs("PNG", output + list[i]);
        close();
    }
}

print("Done!");
