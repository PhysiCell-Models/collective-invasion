import os
import subprocess
import shutil

# assumes that there is a movie making script in each directory - and that the output directory is created
# make sure the movie script points to the python_imaging directory - at the corect level

def run_secondary_scripts(root_dir, output_dir):
    for root, dirs, files in os.walk(root_dir):
        # print(root)
        # print(dirs)
        # print(files)
        for file in files:
            # print(file)
            # if file.endswith(".py") and file != "run_movie_scripts_orthogonal.py" || file != "data_analysis_Painter.py" || file != "pyMCDS_ECM.py"):
            if file == "simple_test_movies_cells_and_environment_Painter.py":
                print(file)
                script_path = os.path.join(root, file)
                folder_name = os.path.basename(root)
                output_mp4 = "orthogonal_" + folder_name + ".mp4"
                print(script_path)
                print(folder_name)
                print(output_mp4)

        #         # Run the secondary script
                os.chdir(folder_name)
                # print(os.getcwd()) 
                
                    # run the script
                subprocess.call(["python", file])
                mp4_path = os.path.join(root, output_mp4)
                # print(mp4_path)

                    # Could confirm is the movie was made - you would start with something like this
                # if os.path.exists(mp4_path):
                    # Rename the .mp4 file to match the folder name
                # new_mp4_path = os.path.join(root, folder_name + ".mp4")
                
                    # give the movie a new name
                os.rename('multi_color_movie.mp4', output_mp4)
                print(f"Renamed multi_color_movie.mp4 to {output_mp4}")
                    
                    # Copy the .mp4 file to the output directory
                output_mp4_copy = os.path.join(output_dir, "orthogonal_" + folder_name + ".mp4")
                shutil.copy(output_mp4, output_mp4_copy)
                print(f"Copied {output_mp4} to {output_mp4_copy}")
                    # figure this one out next. Ideally copy the file in also ... last step 

                os.chdir('../')
                # print(os.getcwd())

if __name__ == "__main__":
    root_directory = "." 
    output_directory = "../../../1_images_for_paper/stochastic_replicates/invasive_cellular_front_orthogonal/" # remember - the script is actually workign on directory lower
    run_secondary_scripts(root_directory, output_directory)

        # could make a more general script that handled all the scripts or do something that copied them all (or just one) into the 
        # directory and then ran them.
    # script_names = ["script1.py", "script2.py", "script3.py"]
    # subdirectory_prefix = "common_prefix"  # Replace with your common prefix
    # for script_name in script_names:
    #     copy_and_run_secondary_scripts(root_directory, script_name, subdirectory_prefix, output_directory)
