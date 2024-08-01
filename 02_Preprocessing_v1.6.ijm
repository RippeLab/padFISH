// This script uses maximum projected tif files and their metadata file and performs flatfield and chromatic aberration correction on them, as well as stitching
// The user will be asked if stitching should be performed

// REQUIRED: input filenames have to be called "MAX_ ... .tif" and have to be in the same folder as the metadata file!

// Collect all choices that user has to make in this script
// Ask user for the data format: single images, images from position list, tilescan
#@ String (choices={"Tilescan", "Multiposition"}, style="radioButtonHorizontal", label="Choose the format of your imaging data") formatchoice

// Ask user if all channels should be used for stitching
#@ Boolean (persist=false, value=true, label="Use all channels for stitching if tilescan") stitchchannelchoice
// (if unticked, select which channels to use for each selected directory individually)

// Ask user if all channels detected should be used for flatfield and chromatic aberration correction or if user wants to select individually for each previously selected directory
#@ Boolean (persist=false, value=true, label="Correct flatfield and chromatic aberration in all channels") correctionchoice
// (if unticked, select which channels to use for each selected directory individually)

// Before starting, close everything that is open
close("*");

// Set options
setOption("ExpandableArrays", true);
setBatchMode(true);

// Define variables
dirnames = newArray;
outdirnames = newArray;
logfilenames = newArray;
logfilelocations = newArray;
var n = 0;

// Generate timestamp (e.g. for output directory name, ...)
function gettimestamp() {
	getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
	ShortTimeString = "" + year;
	month = month + 1;
	if (month < 10) {ShortTimeString = ShortTimeString + "0";}
	ShortTimeString = ShortTimeString + month;
	if (dayOfMonth<10) {ShortTimeString = ShortTimeString+"0";}
	ShortTimeString = ShortTimeString+dayOfMonth;
	if (hour<10) {ShortTimeString = ShortTimeString+"0";}
	ShortTimeString = ShortTimeString+hour;
	if (minute<10) {ShortTimeString = ShortTimeString+"0";}
	ShortTimeString = ShortTimeString+minute;
	if (second<10) {ShortTimeString = ShortTimeString+"0";}
	ShortTimeString = ShortTimeString+second;
	return ShortTimeString;
}

// Generate more detailed timestamp (e.g. for documentation in log file, ...)
function getverbosetimestamp() {
	MonthNames = newArray("01","02","03","04","05","06","07","08","09","10","11","12");
	DayNames = newArray("Sun", "Mon","Tue","Wed","Thu","Fri","Sat");
	getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
	TimeString = DayNames[dayOfWeek]+" ";
	if (dayOfMonth<10) {TimeString = TimeString+"0";}
	TimeString = " "+TimeString+dayOfMonth+"-"+MonthNames[month]+"-"+year+"  ";
	if (hour<10) {TimeString = TimeString+"0";}
	TimeString = TimeString+hour+":";
	if (minute<10) {TimeString = TimeString+"0";}
	TimeString = TimeString+minute+":";
	if (second<10) {TimeString = TimeString+"0";}
	TimeString = TimeString+second;
	return TimeString;
}

timestamp = gettimestamp();
verbosetimestamp = getverbosetimestamp();

// Define function to . . .
function SetLUTResetContrast(LUT) {
	// selectWindow(stack);
	Stack.getDimensions(width, height, channels, slices, frames);
	for(j = 0; j < channels; j++) {
		Stack.setChannel(j+1);
		run(LUT);
		setMinAndMax(0, 65535);				
	}
}

// Select and get directories that should be processed (+ create an output directory and logfile in each selected directory)
while (true) {
	dir = getDirectory("Choose a directory containing your .tif files and their metadata file");
	if (dir == "")
        break; // Exit the loop if the user cancels the directory selection
    dirnames[n] = dir;
    
    // Create output directory
    File.makeDirectory(dirnames[n] + "_" + timestamp + "_preproccessed/");
    outdir = dirnames[n] + "_" + timestamp + "_preproccessed/";
    outdirnames[n] = outdir;
    
    // Create log file
    logfile = File.open(outdirnames[n] + "log.txt");
    logfilenames[n] = logfile;
    logfilelocations[n] = outdirnames[n] + "log.txt";
    print(logfilenames[n], "Run started:       " + verbosetimestamp);
    print(logfilenames[n], "Analyzing data in:  " + dirnames[n]);
    print(logfilenames[n], "");
    File.close(logfilenames[n]);
    
    // Ask if the user wants to process more folders
    test = getString("Add more folders? (y/n)", "y");
    if (test == "n")
    	break; // Exit the loop if the user chooses not to add more folders  
    n = n + 1;
}

// Check directories in the log window
print("Following directories selected for processing:");
for (i = 0; i < dirnames.length; i++) {
    print("  " + dirnames[i]);
}

// Loop over all selected folders and get filelist for each
for (i = 0; i < dirnames.length; i++) {
	filelist = getFileList(dirnames[i]);
	print("Currently processing: " + dirnames[i]);
	
	// Define variables for this run of the loop
	meta = newArray;
	laser_on = newArray;
	wavelengths = newArray;
	correct = newArray;
	ffgain = newArray;
	
	// Check that there is a metafile in the directory and ask the user to select one in case there are multiple (+ add to log file which file is used)
	for (j = 0; j < filelist.length; j++) {
		if (endsWith(filelist[j], "metadata.txt")) {
			meta = Array.concat(meta,filelist[j]);
		}
	}
	if (meta.length == 0) {
		exit("No metadata file found.");
	} else if (meta.length > 1) {
		Dialog.create("Multiple metadata files were found. Select one to extract the channel, objective and montage information from.");
		Dialog.addChoice("metadata file:", meta);
		Dialog.show();
		meta = Dialog.getChoice();
		// Add to log file
		File.append("- Using " + meta + " to extract channel information.", logfilelocations[i]);
		meta = dirnames[i] + meta;
	} else if (meta.length == 1) {
		// Add to log file
		File.append("- Using " + meta[0] + " to extract channel information.", logfilelocations[i]);
		meta = dirnames[i] + meta[0];
	}
	
	// Extract information from the metadata file (channel and stitching configuration, objective, ...)
	ctrl = 0;
	montage = 0;
	objective = 0;
	camera = 0;
	lns = split(File.openAsString(meta), "\n");
	
	for (j = 0; j < lns.length; j++){
		lns[j] = replace(lns[j], "\t", "");
		
		// Extract existing laser lines and their wavelengths
		if (startsWith(lns[j], "{DisplayName=Laser Wavelength")) {
			lambda = split(lns[j],"=");
			laser = replace(lambda[1], "Laser Wavelength ", "");
			laser = replace(laser, ", Value", "");
			laser = parseInt(laser);
			lambda = replace(lambda[2], "}", "");
			lambda = parseInt(lambda);
			List.set(laser, lambda);
		} else if (lns[j] == "[Imaging Modes In Protocol]") {
			ctrl = 1;
		} else if (lns[j] == "[FieldMontageProtocolSpecification]") {
			montage=1;
		}
		
		// Extract used laser lines and stitch parameters (i.e. tiling dimensions)
		if (startsWith(lns[j], "{DisplayName=Laser") && endsWith(lns[j], "True}") && ctrl == 0) {
			temp = replace(lns[j], "\\{DisplayName=Laser ", "");
			temp = replace(temp, ", Value=True}", "");
			laser_on = Array.concat(laser_on, parseInt(temp));
		} else if (startsWith(lns[j], "Rows") && montage == 1) {
			temp = replace(lns[j],"Rows=", "");
			y_grid = temp;
		} else if (startsWith(lns[j], "Columns") && montage == 1) {
			temp = replace(lns[j], "Columns=", "");
			x_grid = temp;
		} else if (startsWith(lns[j], "Overlap") && montage == 1) {
			temp = replace(lns[j], "Overlap=", "");
			o = temp;
		}
		
		// Extract the objective information from the metadata file
		if (startsWith(lns[j], "{DisplayName=Select Objective") && objective == 0) {
			temp = split(lns[j], " - ");
			objective = split(temp[4], "/");
			objective = objective[0];		
		}
		
		// Extract the camera information from the metadata file
		if (startsWith(lns[j], "AcquisitionDeviceAlias=") && camera == 0) {
			temp = split(lns[j], " ");
			camera = temp[1] + temp[2];
		}
	}
	
	//remove first element of laser_on --> when using the twin cam, the Zyla laser is listed twice as on
	//laser_on = Array.slice(laser_on,1);
	
	// Get list with all existing lasers and their corresponding wavelengths (+ add to log file)
	laser = List.getList;
	for (j = 0; j < laser_on.length; j++) {
		wavelengths = Array.concat(wavelengths,List.get(laser_on[j]));
	}
	// Add to log file
	//File.append("- Found " + wavelengths.length + " active channels with the following excitation wavelengths:", logfilelocations[i]);	
	//for (j = 0; j < wavelengths.length; j++) {
		// Add to log file
	//	File.append("  " + wavelengths[j] + " nm", logfilelocations[i]);
	//}
	
	// Stitch if selected
	if (formatchoice == "Tilescan") {
		if (stitchchannelchoice == 1) {
			// Automatically select all channels for stitching
			Dialog.create("Select channel(s) to be used for stitching.\n If multiple channels are selected, a maximum intensity projection of them will be used.");
			for (j = 0; j < wavelengths.length; j++) {
				Dialog.addCheckbox(wavelengths[j] + " nm", true);
			}
			for (j = 0; j < wavelengths.length; j++) {
				stitch = Array.concat(stitch,Dialog.getCheckbox());
			}
		} else {
			// Choose the channel to be used for stitching (usually DAPI)
			Dialog.create("Select channel(s) to be used for stitching.\n If multiple channels are selected, a maximum intensity projection of them will be used.");
			for (j = 0; j < wavelengths.length; j++) {
				Dialog.addCheckbox(wavelengths[j] + " nm", true);
			}
			Dialog.show();
			for (j = 0; j < wavelengths.length; j++) {
				stitch = Array.concat(stitch,Dialog.getCheckbox());
			}
		}
		// Remove additional array value (first)
		stitch = Array.deleteIndex(stitch, 0);
		
		// Get the index and the wavelengths of the channel(s) to be used for stitching
		for (j = 0; j < wavelengths.length; j++) {
			if (stitch[j] == 1) {
				stitch_index = Array.concat(stitch_index, j+1);
				stitch_channel = Array.concat(stitch_channel, wavelengths[j]);
			}
		}
		
		stitch_index = Array.deleteIndex(stitch_index, 0);
		stitch_channel = Array.deleteIndex(stitch_channel, 0);
		stitch_print = String.join(stitch_channel, ",");
		// Add to log file
		File.append("- Stitching a "+ x_grid + " x " + y_grid + " image with " + o + "% overlap between the tiles based on the " + stitch_print + " nm channel(s).", logfilelocations[i]);
	} else {
		// Add to log file
		File.append("- No stitching was performed.", logfilelocations[i]);
	}
	
	if (correctionchoice == 1) {
		// Automatically select all channels to be corrected (flatfield AND chromatic aberration)
		Dialog.create("Flatfield and chromatic aberration correction:");
		for (j = 0; j < wavelengths.length; j++) {
			Dialog.addCheckbox(wavelengths[j] + " nm", true);
		}
		for (j = 0; j < wavelengths.length; j++) {
			correct = Array.concat(correct,Dialog.getCheckbox());
		}
	} else {
		// Select channels that should be corrected (flatfield AND chromatic aberration)
		Dialog.create("Flatfield and chromatic aberration correction:");
		for (j = 0; j < wavelengths.length; j++) {
			Dialog.addCheckbox(wavelengths[j] + " nm", true);
		}
		Dialog.show();
		for (j = 0; j < wavelengths.length; j++) {
			correct = Array.concat(correct,Dialog.getCheckbox());
		}
	}
	
	// Add to log file
	//File.append("Correcting the flatfield using the " + objective + " correction files for " + camera, logfilelocations[i]);
	
	// Get gain files for flatfield correction
	correctdir = "/Volumes/sd17B002/Microscopy-Dragonfly/Image-Analysis/Flat-Field-Correction/2022-05-16/";
	
	for (j = 0; j < wavelengths.length; j++) {
		gainfile = correctdir + camera + "-" + wavelengths[j] + "nm-" + objective + "-2022-05-16_GainImage.tif";
		ffgain = Array.concat(ffgain, gainfile);
	}
	ffdark = correctdir + camera + "darkCounts.tif";
	
	// Get transform folder for chromatic aberration correction
	transformdir = "/Volumes/sd17B002/Microscopy-Dragonfly/Image-Analysis/ChromaticAbberationCorrection/transform-XML-file-to405nm-beads1um_27Nov2023/";
	
	// Add to log file
	File.append("- Using the following files for flatfield correction:", logfilelocations[i]);
	for (j = 0; j < ffgain.length; j++) {
		File.append("  " + ffgain[j], logfilelocations[i]);
	}
	File.append("  " + ffdark, logfilelocations[i]);
	File.append("- Using the following folder for chromatic aberration correction:", logfilelocations[i]);
	File.append("  " + transformdir, logfilelocations[i]);
	
	// Create further output directories
    tempdir = outdirnames[i] + "_temp/";
    outdirmax = outdirnames[i] + "_TIF_MaxProj/"; //"_max_tif/";
    fused = outdirnames[i] + "_Fused/";
    stitchdir = outdirmax + "_Stitch/";
    if (formatchoice == "Tilescan") {
		File.makeDirectory(fused);
	}
    File.makeDirectory(tempdir);
    File.makeDirectory(tempdir + "out/");
	File.makeDirectory(outdirmax);
	File.makeDirectory(stitchdir);
	
	for (j = 0; j < wavelengths.length; j++) {
		File.makeDirectory(outdirmax + "/_" + wavelengths[j] + "/");
		if (correct[j] == 1) {
			File.makeDirectory(outdirmax + "/_" + wavelengths[j] + "/_corrected");
		}
	}
	
	// Empty stitch folder (in case macro was run before)
	if (formatchoice == "Tilescan") {
		templist = getFileList(stitchdir);
		for (j = 0; j < templist.length; j++) {
			File.delete(stitchdir + templist[j]);
		}
	}
	
	// Count number of .ims files in directory
    nimsfiles = 0;
    for (j = 0; j < filelist.length; j++) {
    	if (endsWith(filelist[j], ".ims")) {
    		nimsfiles = nimsfiles + 1;
        }
    }
    o = 1;
	
	// Perform maximum projection, flatfield and chromatic aberration correction
	for(j = 0; j < filelist.length; j++) {
		if (endsWith(filelist[j], ".ims")) {
			filenm = dirnames[i] + filelist[j];
			check = dirnames[i] + "MAX_" + replace(filelist[j],".ims",".tif");
			tifname = "MAX_" + replace(filelist[j],".ims",".tif");
			//print(check);
			if (File.exists(check)) {
				open(check);
				rename("MAX_Stack");
				print("Tif available for .ims file " + o + "/" + nimsfiles + " (" + tifname + ")");
				File.append("  Tif available for .ims file " + o + "/" + nimsfiles + " (" + tifname + ")", logfilelocations[i]);
			} else {
				print("Converting .ims file " + o + "/" + nimsfiles + " (" + filelist[j] + ")");
				File.append("  Converting .ims file " + o + "/" + nimsfiles + " (" + filelist[j] + ")", logfilelocations[i]);
				run("Bio-Formats (Windowless)", "open=filenm");
				SetLUTResetContrast("Grays");
				rename("3DStack");
				run("Z Project...", "projection=[Max Intensity]");
				setMinAndMax(0, 65535);
				saveAs("Tiff", dirnames[i] + "MAX_" + filelist[j]);
				rename("MAX_Stack");
				close("3DStack");
			}
			for (k = 0; k < wavelengths.length; k++) {
				channel_index = k+1;
				outdirmaxch = outdirmax + "/_" + wavelengths[k] + "/";
				selectImage("MAX_Stack");
				run("Duplicate...", "duplicate channels=channel_index-channel_index");
				setMinAndMax(0, 65535);
				saveAs("Tiff", outdirmaxch + wavelengths[k] + "nm_MAX_" + filelist[j]);
				rename("Raw");
				if (correct[k] == 0) {
					print("No correction for " + wavelengths[k] + " nm images. Continuing with next channel...");
					close("Raw");
					continue;
				}
				print("    Correcting flatfield and chromatic aberration for " + wavelengths[k] + " nm channel of " + tifname);
				File.append("    Correcting flatfield and chromatic aberration for " + wavelengths[k] + " nm channel of " + tifname, logfilelocations[i]);
				
				// Calculate the corrected Image (RawImage-Dark) x Gain
				open(ffgain[k]);
				rename("gain");
				setMinAndMax(0, 65535);
				open(ffdark);
				rename("dark");
				setMinAndMax(0, 65535);
				imageCalculator("Subtract create", "Raw", "dark");
				rename("raw_minus_dark");
				run("32-bit");
				setMinAndMax(0, 65535);
				imageCalculator("Multiply create 32-bit", "raw_minus_dark" , "gain");
				rename("corr_Image");  //flatfield-corrected image
				setMinAndMax(0, 65535);
				run("16-bit");
				setMinAndMax(0, 65535);
				tempimage = wavelengths[k] + "nm_MAX_corr_" + replace(filelist[j],".ims",".tif");
				saveAs("Tiff", tempdir + tempimage);
				close();
				close("Raw");
				close("gain");
				close("dark");
				close("raw_minus_dark");
				run("Transform Virtual Stack Slices", "source=" + tempdir + " output=" + tempdir + "out transforms=" + transformdir + camera + "_" + wavelengths[k] + "nm_" + objective + " interpolate");
				close();
				open(tempdir + "out/" + tempimage);
				if (wavelengths[k] == 730 && camera == "EMCCD1" && objective == "100x") {
					makeRectangle(0, 1, 1024, 1024);
					run("Crop");
				}
				saveAs("Tiff", outdirmaxch + "_corrected/" + tempimage);
				close();
				discard = File.delete(tempdir + tempimage);
				discard = File.delete(tempdir + "out/" + tempimage);
			}
			close("*");
			
			// Continue with stitching if selected by user
			if (formatchoice == "Tilescan") {
				// Open all corrected images to be used for stitching
				for (k = 0; k < stitch_channel.length; k++) {
					open(outdirmax + "_" + stitch_channel[k] + "/_corrected/" + stitch_channel[k] + "nm_MAX_corr_" + replace(filelist[j],".ims",".tif"));
					//run("Enhance Contrast", "saturated=0.01");
					//run("Apply LUT");
					setMinAndMax(0, 65535);
				}
				if (stitch_channel.length > 1) {
					run("Images to Stack", "use");
					run("Z Project...", "projection=[Max Intensity]");
				}
				saveAs("Tiff", stitchdir + String.join(stitch_channel, "_") + "nm_MAX_corr_" + replace(filelist[j], ".ims", ".tif"));
				close("*");
			}
			o = o + 1;
		} else {
			//print("Skipped "+filelist[j]);
		}
	}
	File.append("- Finished correction for all files.", logfilelocations[i]);
	
	// Continue with stitching if selected by user
	if (formatchoice == "Tilescan") {
		// Extract index of first image (..._F????.tif)
		maxlist = getFileList(stitchdir);
		filenm = maxlist[0];
		init = split(maxlist[0], "_");
		init = init[init.length-1];
		init = replace(init, "F", "");
		init = replace(init, ".tif", "");
		init = parseInt(init);
		
		filenm = replace(filenm , "[0-9]{4}.tif", "\\{iiii\\}.tif");
		filenm = replace(filenm , "[0-9]{3}.tif", "\\{iii\\}.tif");
		filenm = replace(filenm , "[0-9]{2}.tif", "\\{ii\\}.tif");
		filenm = replace(filenm , "[0-9]{1}.tif", "\\{i\\}.tif");
		
		// Run stitching on selected channel (computing perfect overlap)
		outdirmaxstitch = outdirmax + "/_stitch/";
		run("Grid/Collection stitching", "type=[Grid: snake by columns] order=[Up & Right] grid_size_x=x_grid grid_size_y=y_grid tile_overlap=o first_file_index_i=init directory=" + outdirmaxstitch + " file_names=" + filenm + " output_textfile_name=TileConfiguration.txt fusion_method=[Do not fuse images (only write TileConfiguration)] regression_threshold=0.30 max/avg_displacement_threshold=2.5 absolute_displacement_threshold=3.5 compute_overlap subpixel_accuracy computation_parameters=[Save computation time (but use more RAM)] image_output=[Fuse and display]");
		
		// Modify file names in the output tile configuration file for each channel
		lns = split(File.openAsString(outdirmaxstitch + "TileConfiguration.registered.txt"), "\n");
		
		// in uncorrected maximum projection folder (stack with all colors)
		outfile = outdirmax + "TileConfiguration.registered.txt";
		if (File.exists(outfile)) {
			File.delete(outfile);
		}
		for (j = 0; j < lns.length; j++) {
			if (startsWith(lns[j], String.join(stitch_channel, "_"))) {
				lns[j] = replace(lns[j], String.join(stitch_channel, "_") + "nm_MAX_corr", "MAX");
				File.append(lns[j], outfile);
			} else {
				File.append(lns[j], outfile);
			}
		}
		for (j = 0; j < wavelengths.length; j++) {
			//in uncorrected maximum projection folder (each color)
			outdirmaxch = outdirmax + "_" + wavelengths[j] + "/";
			outfile2 = outdirmaxch + "TileConfiguration.registered.txt";
			if (File.exists(outfile2)) {
				File.delete(outfile2);
			}
			temp_lns = Array.copy(lns);
			for (k = 0; k < temp_lns.length; k++) {
				if (startsWith(temp_lns[k], "MAX")) {
					temp_lns[k] = replace(temp_lns[k], "MAX", wavelengths[j] + "nm_MAX");
					File.append(temp_lns[k], outfile2);
				} else {
					File.append(temp_lns[k], outfile2);
				}
			}
			
			// in corrected maximum projection folder (if correction was selected)
			if (correct[j] == 0) {
				continue;
			}
			outdirmaxchcorr = outdirmaxch + "_corrected/";
			outfile3 = outdirmaxchcorr + "TileConfiguration.registered.txt";
			if (File.exists(outfile3)) {
				File.delete(outfile3);
			}
			temp_lns = Array.copy(lns);
			for (k = 0; k < temp_lns.length; k++) {
				if (startsWith(temp_lns[k], "MAX")) {
					temp_lns[k] = replace(temp_lns[k], "MAX", wavelengths[j] + "nm_MAX_corr");
					File.append(temp_lns[k], outfile3);
				} else {
					File.append(temp_lns[k], outfile3);
				}
			}
		}
		
		// Apply stitch positions to uncorrected and (if applicable) corrected single-color images
		filenm = replace(filenm, String.join(stitch_channel, "_") + "nm_MAX_corr", "MAX");
		filenm = replace(filenm, "_F\\{i{1,4}\\}.tif", "_Fused_");
		for (j = 0; j < wavelengths.length; j++) {
			outdirmaxch = outdirmax + "_" + wavelengths[j] + "/";
			outdirmaxchcorr = outdirmaxch + "_corrected/";
			run("Grid/Collection stitching", "type=[Positions from file] order=[Defined by TileConfiguration] directory=" + outdirmaxch + " layout_file=TileConfiguration.registered.txt fusion_method=[Linear Blending] regression_threshold=0.30 max/avg_displacement_threshold=2.50 absolute_displacement_threshold=3.50 subpixel_accuracy computation_parameters=[Save computation time (but use more RAM)] image_output=[Fuse and display]");
			setMinAndMax(0, 65535);
			saveAs("Tiff", fused+filenm+wavelengths[j]+"nm.tif");
			if (correct[j] == 0) {
				continue;
			}
			run("Grid/Collection stitching", "type=[Positions from file] order=[Defined by TileConfiguration] directory=" + outdirmaxchcorr + " layout_file=TileConfiguration.registered.txt fusion_method=[Linear Blending] regression_threshold=0.30 max/avg_displacement_threshold=2.50 absolute_displacement_threshold=3.50 subpixel_accuracy computation_parameters=[Save computation time (but use more RAM)] image_output=[Fuse and display]");
			setMinAndMax(0, 65535);
			saveAs("Tiff", fused+filenm+wavelengths[j]+"nm_corr.tif");
		}
		close("*");		
	}
	discard = File.delete(tempdir + "out/");
	discard = File.delete(tempdir);
	
	// Delete stitch directory if not used
	if (formatchoice == "Multiposition") {
		discard = File.delete(stitchdir);
	}
	// Add to log file
	//File.append("- Finished " + outdirnames[i], logfilelocations[i]);
	File.append("", logfilelocations[i]);
	File.append("Run finished:      " + getverbosetimestamp(), logfilelocations[i]);
}

// Finish it off
if (formatchoice == "Tilescan") {
	print("Finished - completed corrections and stitching");
} else {
	print("Finished - completed corrections");
}
close("*");
run("Collect Garbage");