# Learning rate per variant: the one with the lower validation loss at seed 0, width 128. Writes results/settings.json.
import json
from pathlib import Path

import pandas as pd

rows = []
for path in Path("results").glob("*/seed0_width128_lr*/run.json"):
    run = json.load(open(path))
    rows.append({"variant": path.parent.parent.name, "mask": run["mask"], "lr": run["lr"],
                 "val_loss": round(run["val_loss"], 1), "best_epoch": run["best_epoch"]})
table = pd.DataFrame(rows).sort_values(["variant", "lr"])
print(table.to_string(index=False))
assert (table.groupby("variant").size() == 2).all(), "every variant needs a run at both learning rates"

chosen = table.loc[table.groupby("variant")["val_loss"].idxmin()]
settings = {row.variant: {"mask": row.mask, "lr": row.lr} for row in chosen.itertuples()}
print(json.dumps(settings, indent=2))
with open("results/settings.json", "w") as f:
    json.dump(settings, f, indent=2)
