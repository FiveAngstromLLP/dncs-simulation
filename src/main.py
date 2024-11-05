
import os
import toml
import time
import dncs
import shutil
from typing import List
from dataclasses import dataclass
from integrator import DncsIntegrator, CleanUp, MDSimulation


@dataclass
class SimulationConfig:
    sequence: str
    moleculename: str
    interface: str
    n_samples: int
    temp: float
    forcefield: List[str]
    device: str
    solvent: int
    steps: int
    gamma: float
    dt: float
    md_steps: int
    grid: int

class GenerateSamples:
    def __init__(self, config: SimulationConfig):
        self.config = config
        self.polymer = None
        self.sample = None

    def generate_samples(self):
        self.polymer = dncs.Polymer(self.config.sequence, "amberfb15.xml")
        self.sample = dncs.SobolSampler(self.polymer, self.config.n_samples, self.config.grid)


    def save_to_pdb(self):
        folder = os.environ.get('DNCS_FOLDER')
        if os.path.exists(f"{folder}/Result/{self.config.moleculename}"):
            shutil.rmtree(f"{folder}/Result/{self.config.moleculename}")
        os.makedirs(f"{folder}/Result/{self.config.moleculename}")
        self.sample.toPDB(f"{folder}/Result/{self.config.moleculename}/multi_structure.pdb")
        self.sample.toPDBFiles(f"{folder}/Result/{self.config.moleculename}/Sampled")
        self.sample.write_angles(f"{folder}/Result/{self.config.moleculename}/Sampled/angles.out")

def runopenmm(config):
    # Run Integrator
    integrator = DncsIntegrator(config)
    integrator.run_integrator()
    CleanUp(config)

    md = MDSimulation(config)
    md.run_simulation()



if __name__ == "__main__":
    # Load configuration from dncs.toml file
    with open('dncs.toml', 'r') as f:
        toml_config = toml.load(f)

    config = SimulationConfig(**toml_config['simulation'])

    start_time = time.time()


    # Generate conformational samples
    conformation = GenerateSamples(config)
    conformation.generate_samples()
    conformation.save_to_pdb()

    if config.interface == "openmm":
        runopenmm(config)

    end_time = time.time()
    simulation_time = end_time - start_time

    if simulation_time > 60:
        minutes = simulation_time / 60
        print(f"Total simulation time: {minutes:.5f} minutes")
    else:
        print(f"Total simulation time: {simulation_time:.5f} seconds")

    # Final cleanup
    import gc
    gc.collect()
