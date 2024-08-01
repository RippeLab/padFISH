// this is a macro that can run in both headless and GUI mode
// in headless mode it uses parameters defined from command line/a shell script 

#@ File (style="directory", persist=false, label="Image directory") dir
#@ Boolean (persist=false, label="Save z-stack (unprojected)?") saveZstack
#@ Boolean (persist=false, label="Save Maximum projection?") saveMaxproj
#@ Boolean (persist=false, value=false, label="Save cropped image for dSTORM?") saveCrop
#@ String (persist=false, label="Image resolution", value="1024x1024", choices={"1024x1024", "2048x2048"}, style="radioButtonHorizontal") imres

dir = dir+"/";

// get directories from a tool box
list = getFileList(dir);

setBatchMode(true);
//do the tiff-converstation for ims files in the directory
for(i=0; i<list.length; i++) {
	if (endsWith(list[i], ".ims")){
		filenm = dir+list[i];
		print("Processing file "+(i+1)+"/"+list.length);
		run("Bio-Formats Macro Extensions");
		Ext.openImagePlus(filenm);
		
		if (saveZstack == true){
			saveAs("Tiff", dir+list[i]);
			print("saved Z-stack .tif of " +list[i]);
		}
		
		if (saveMaxproj == true){
			run("Z Project...", "projection=[Max Intensity]");
			saveAs("Tiff", dir+ "MAX_" +list[i]);
			print("saved Maxproj .tif of " +list[i]);
		}
		
		if (saveCrop == true){
			
			if (imres == "1024x1024"){
				makeRectangle(330, 410, 325, 325);
			} else if (imres == "2048x2048"){
				makeRectangle(660, 820, 650, 650);
			}
			run("Crop");
			saveAs("Tiff", dir+ "Cropped_" +list[i]);
			print("saved cropped .tif of " +list[i]);
		}
		
		close("*");
	}else{
		print("skipped "+list[i]);
	}
}
print("completed conversion");
