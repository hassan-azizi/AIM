# MATLAB Project Template File for Compiling Source Code with MATLAB

This directory contains a MATLAB project file (`AIM_Build.prj`) designed to facilitate end users who want to compile the provided source code using the MATLAB for re-distribution on machines with OS other than Microsoft Windows. The successful compilation will generate AIM installer that can be distributed.
Follow the instructions below to set up and compile the project.

## Prerequisites
- **MATLAB Version**: Ensure you have MATLAB R2024a or a newer version installed.
- **MATLAB Compiler**: This MATLAB toolbox is required to compile MATLAB code into standalone applications. Version must be R2024a or newer.
- **MATLAB Runtime**: The compiled application requires the MATLAB Runtime environment to run on machines without MATLAB installed. The version of MATLAB Runtime must match the version of MATLAB Compiler used for compilation. While compiling the source code, one has two options
    1. Create an application installer with MATLAB Runtime.
    2. Create an application installer without MATLAB Runtime. In this case, users of the application will need to download the MATLAB Runtime.
<br>In the former case, you must download the MATLAB Runtime installer that matches both the version and update level of MATLAB used to create the installer.<br>
Details about installing MATLAB Runtime on different platforms can be found [here](https://www.mathworks.com/help/compiler/install-the-matlab-runtime.html).

## Setup Instructions

1. **Clone or Download the Repository**:
   - Clone the repository using Git:
     ```bash
     git clone https://github.com/mtap-research/AIM.git
     ```
   - Or download the ZIP file and extract it to your desired location.

2. **Start MATLAB**:
   - Start the MATLAB application.

3. **Open the AIM Project**:
   - Navigate to the AIM/build directory.  
   - Double-click the `AIM_Build.prj` file to open it in MATLAB.

4. **Add Necessary Files**:
   - In the MATLAB Compiler project window, add all the files present in the AIM/src directory.

5. **Configure Application Settings**:
   - Specify the main file to be executed, by adding "Main_app.mlapp" from AIM/src directory.

6. **Compile the Application**:
   - Click the "Package" button to start the compilation process.
   - Upon completion, MATLAB will generate the standalone application and a `readme.txt` file containing deployment prerequisites and a list of packaged files.

## Deployment

- **MATLAB Runtime Installation**:
  - Users must install the MATLAB Runtime that matches the version of MATLAB Compiler used. The installer can be downloaded from the MathWorks website.

- **Running the Application**:
  - After installing MATLAB Runtime, users can run the compiled application by executing the generated executable file.

## Additional Resources

- **MATLAB Compiler Documentation**: Detailed information on using the MATLAB Compiler can be found [here](https://www.mathworks.com/help/compiler/getting-started-with-matlab-compiler.html).
- **MATLAB Runtime Documentation**: Learn more about MATLAB Runtime [here](https://www.mathworks.com/products/compiler/matlab-runtime.html).

For any issues or questions, please contact us at [hassanaz.14@outlook.com].
