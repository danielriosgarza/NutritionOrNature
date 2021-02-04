# NutritionOrNature

Repository of Figures and scripts used in the manuscript: [Nutrition or nature: disentangling the complex forces shaping prokaryote pan-genomes](https://www.biorxiv.org/content/10.1101/2020.12.14.422685v3)

Below is a guide to enable readers to reproduce our results:

1) The environment ball:

We simulated pan-reactomes in a random ball containing the relative concentration of metabolites. The environments we used in the manuscript are listed in [TableS3](https://github.com/danielriosgarza/NutritionOrNature/blob/main/Table_S3.xlsx). We used the script: [generate_env_ball.py](https://github.com/danielriosgarza/NutritionOrNature/blob/main/Scripts/generate_env_ball.py) to generate these environments.

An example of the expected result of this script is available [here](https://github.com/danielriosgarza/NutritionOrNature/blob/main/Files/env_ball_1000.tsv)

The random environment generated is handled by importing the [Env_ball_class.py](https://github.com/danielriosgarza/NutritionOrNature/blob/main/Scripts/Env_ball_class.py) script, for example:

```python
from Env_ball_class import Env_ball
ev = Env_ball(1000)
ev.plot()
print(ev.metabolites)

```

2) Create family-level pan-reactomes:

We begin with a database of draft reconstractions of strain-level genome-scale metabolic models (GSMMs) listed in [TableS2](https://github.com/danielriosgarza/NutritionOrNature/blob/main/Table_S2.xlsx). The models were reconstructed using [ModelSEED](https://modelseed.org/) and are available on request.

To create family-level pan-reactomes we used the script [create_panReactome_model.py](https://github.com/danielriosgarza/NutritionOrNature/blob/main/Scripts/create_panReactome_model.py). 

Usage example:
```python
from create_panReactome_model import PanReactome_model
em = PanReactome_model('Aeromonadaceae',PathToModels)
panReac = em.make_panReactome()
em.add_transporter(panReac)
em.write_model(panReac, outputModel)

```

An example of a model reconstructed for the Aeromonadaceae family that is extensively described in the manuscript is available here: [Aeromonadaceae_reactome](https://github.com/danielriosgarza/NutritionOrNature/blob/main/Files/Aeromonadaceae.ensembl.sbml). 


3) The toy model example

The code for building and simulating the toy model is available through the script [ToyModelSimulator.py](https://github.com/danielriosgarza/NutritionOrNature/blob/main/Scripts/toyModelSimulator.py). The script can be runned directly in a python environment and generates a pickle that contains a dictionary named "store" with all the relevant information about the toy model simulation. The script takes some hours to complete (due to the Moran process).

4) generate panEFMs

To generate panEFMs from a pan-reactome (generated in step2) use the script [get_PanEFMs.py](https://github.com/danielriosgarza/NutritionOrNature/blob/main/Scripts/get_PanEFMs.py). 
This script also takes several hours to complete, and example of usage with 20 preocessor jobs follows:

```python

eb=Env_ball(1000)
transporters=eb.transporters[:]
random_environments=eb.matrix.copy()
model = cobra.io.read_sbml_model(PathToPanReactomeModel)
model.solver='gurobi'
reactions=[i.id for i in model.reactions if 'rxn' in i.id]
end_env=job_state*20
start_env=end_env-20
[get_run(model, reactions, transporters, random_environments[i], 'pathToResults'+str(i)+'.tsv', 1000, job_state*3) for i in range(start_env, end_env,1)]

```

4) Analyze the panEFMs

To analyze panEFMs generated from a family-level pan-reactome we used the script [analyzePanEFMs.py](https://github.com/danielriosgarza/NutritionOrNature/blob/main/Scripts/analyzePanEFMs.py). This script first calls the class [parse_panEFM_class.py](https://github.com/danielriosgarza/NutritionOrNature/blob/main/Scripts/parse_panEFM_class.py) that parses the simulation results and gets all the necessary information from strain-level models, next a python dictionary with all the metrics (distances, scores, etc is stored in a pickle).


5) Get Elastic net prediction
An illustrative example of the code we used to build elastic net models is the script [EN_Predictions.py](https://github.com/danielriosgarza/NutritionOrNature/blob/main/Scripts/EN_Predictions.py).
