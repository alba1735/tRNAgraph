import subprocess
import os
from pathlib import Path

class CommandRunner:
    def __init__(self, mode: str = "local", image: str = "trnagraph:latest"):
        """
        mode: 'local', 'docker', or 'singularity'
        """
        self.mode = mode
        self.image = image

    def run(self, cmd_list: list, cwd: str = os.getcwd()):
        if self.mode == "local":
            return self._run_local(cmd_list, cwd)
        elif self.mode == "docker":
            return self._run_docker(cmd_list, cwd)
        elif self.mode == "singularity":
            return self._run_singularity(cmd_list, cwd)

    def _run_local(self, cmd, cwd):
        # Your standard subprocess call
        print(f"Running: {' '.join(cmd)}")
        return subprocess.run(cmd, cwd=cwd, check=True)

    def _run_docker(self, cmd, cwd):
        # We must mount the current working directory into the container
        cwd_path = Path(cwd).resolve()
        
        docker_cmd = [
            "docker", "run", "--rm",
            "-v", f"{cwd_path}:/data",  # Mount current dir to /data
            "-w", "/data",              # Set working dir inside container
            "--user", f"{os.getuid()}:{os.getgid()}", # Match user permissions
            self.image
        ] + cmd
        
        print(f"Docker container call: {' '.join(docker_cmd)}")
        return subprocess.run(docker_cmd, check=True)

    def _run_singularity(self, cmd, cwd):
        # HPCs usually prefer Apptainer/Singularity
        singularity_cmd = ["singularity", "exec", self.image] + cmd
        return subprocess.run(singularity_cmd, cwd=cwd, check=True)