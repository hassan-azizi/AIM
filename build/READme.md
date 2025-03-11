AIM Build
=========================
# MATLAB Project for Compiling Source Code with MATLAB Runtime Compiler

This repository contains a MATLAB project (`*.prj` file) designed to compile the provided source code using the MATLAB Runtime Compiler. Follow the instructions below to set up and compile the project.

## Prerequisites

- **MATLAB Version**: Ensure you have MATLAB R2024a or a newer version installed.
- **MATLAB Compiler**: This toolbox is required to compile MATLAB code into standalone applications. Version must be R2024a or newer.
- **MATLAB Runtime**: The compiled application requires the MATLAB Runtime environment to run on machines without MATLAB installed. The version of MATLAB Runtime must match the version of MATLAB Compiler used for compilation.

## Setup Instructions

1. **Clone or Download the Repository**:
   - Clone the repository using Git:
     ```bash
     git clone https://github.com/yourusername/your-repository.git
     ```
   - Or download the ZIP file and extract it to your desired location.

2. **Open the MATLAB Project**:
   - Navigate to the project directory.
   - Double-click the `YourProjectName.prj` file to open it in MATLAB.

3. **Add Necessary Files**:
   - In the MATLAB Compiler project window, add your main MATLAB files and any required supporting files.

4. **Configure Application Settings**:
   - Specify the main file to be executed.
   - Define the output settings, such as the name of the compiled application.

5. **Compile the Application**:
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

For any issues or questions, please contact [Your Name] at [your.email@example.com].
