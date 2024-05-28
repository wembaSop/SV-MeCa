# Use the Miniconda base image
FROM continuumio/miniconda3:latest

# Set the working directory to /workspace
WORKDIR /workspace
RUN mkdir input

# Install required packages
RUN conda config --add channels conda-forge && conda config --add channels bioconda
RUN conda install -y mamba=1.5.8

RUN mamba install -y nextflow=23.10.1 singularity=2.4.2

# Copy the Python code from the local directory to the container
#COPY align_env.yml /workspace/

#RUN mamba env create -y -f align_env.yml
# Clone the GitHub repository
RUN apt-get update && \
    apt-get install -y git && \
    git clone https://github.com/wembaSop/SV-MeCa.git

#move to the SV-MeCa repository
WORKDIR /workspace/SV-MeCa

# add the tardis image
COPY tardis_latest.img /workspace/SV-MeCa/bin/

#create the required environments
RUN for file in $(ls conda/*.yml);do mamba env create -y -f $file;done

# Set the entrypoint to a shell and change to the workspace directory
ENTRYPOINT ["/bin/bash", "-c"]
CMD ["cd /workspace && /bin/bash"]
