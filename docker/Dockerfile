# Use the Miniconda base image
FROM continuumio/miniconda3:latest

# Set the working directory to /workspace
WORKDIR /workspace
RUN mkdir input

# Install required packages
RUN conda config --add channels conda-forge && conda config --add channels bioconda
RUN conda install -y mamba=1.5.8

RUN mamba install -y nextflow=23.10.1 singularity=3.8.6

# Copy the Python code from the local directory to the container
#COPY align_env.yml /workspace/

#RUN mamba env create -y -f align_env.yml
# Clone the GitHub repository
RUN apt-get update && \
    apt-get install -y git && \
    git clone https://github.com/wembaSop/SV-MeCa.git

#move to the SV-MeCa repository
WORKDIR /workspace/SV-MeCa

# move the tardis image
RUN mv /workspace/SV-MeCa/data/tardis_latest.img /workspace/SV-MeCa/bin/

#create the required environments
RUN for file in $(ls conda/*.yml);do mamba env create -y -f $file;done

# fix version bug of samtools
RUN mamba install -n breakdancer samtools=1.17
RUN mamba install -n pindel samtools=1.17

# Set the entrypoint to a shell and change to the workspace directory
ENTRYPOINT ["/bin/bash", "-c"]
CMD ["cd /workspace && /bin/bash"]
