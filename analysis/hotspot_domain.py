# Mapping hidden urban exposure around Bangladesh's protected areas and species habitats
# Dr. Rakibul Ahasan

import os, sys
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path[:0] = [REPO, os.path.join(REPO, "figures"), os.path.join(REPO, "sdm")]
from config import HARM, MR_TAB, TAB, OUT_FIG, SDM_RASTERS, SDM_WORK, SDM_DATA, REV, SUPP
os.chdir(SDM_WORK)
import numpy as np, pandas as pd
from gistar import gistar
z=np.load("urban_sdm_grid.npz"); L=z["L"]; CA=z["CA"]; irr=np.load("irr.npy"); pri=np.load("pri.npy"); nval=np.load("nval.npy"); rt=np.load("rich_t.npy").astype(float)
iv=irr>-1e30; valid=nval==1098; land=L>0
A=valid&iv&(L/CA>0.5)
rem=iv&~valid&land; ext=land&~iv; thin=land&valid&iv&~A
P=["GLC30","GHSL100","ESA300"]
rows=[]
for nm,m in [("national land",land),("analysis domain",A),("removed as built-up by the models",rem),("outside model extent",ext),("cells less than half land",thin)]:
    r=dict(zone=nm,land_km2=L[m].sum())
    for p in P:
        r[f"{p}_urban2020_km2"]=z[f"{p}_2020"][m].sum(); r[f"{p}_gain_km2"]=(z[f"{p}_2020"]-z[f"{p}_2000"])[m].sum()
    rows.append(r)
acc=pd.DataFrame(rows); acc.to_csv(os.path.join(REV, "domain_accounting.csv"),index=False); print(acc.round(0).to_string(index=False))
print("priority areas: % of domain land", round(100*L[A&(pri==1)].sum()/L[A].sum(),1), "| % of national land", round(100*L[land&(pri==1)].sum()/L[land].sum(),1), "| % of irr-valid cell area", round(100*CA[iv&(pri==1)].sum()/CA[iv].sum(),1))
out=[]
for ver,D in [("analysis domain",A),("all land",land&(L/CA>0.5))]:
    irr_d=np.where(A,irr,0.0); rt_d=np.where(A,rt,0.0)
    RH=gistar(rt_d,D)>2.58; IH=D&(irr_d>=1)
    for k in (3,5):
        RHk=gistar(rt_d,D,k)>2.58
        for p in P:
            UH=gistar(z[f"{p}_2020"]-z[f"{p}_2000"],D,k)>2.58
            for surf,S in [("threatened richness",RHk),("irreplaceability>=1",IH)]:
                obs=L[UH&S].sum()/L[UH].sum(); exp=L[S&D].sum()/L[D].sum()
                out.append(dict(domain=ver,neighborhood=f"{k} by {k}",product=p,surface=surf,domain_km2=L[D].sum(),hotspot_km2=L[UH].sum(),observed_pct=100*obs,expected_pct=100*exp,ratio=obs/exp))
o=pd.DataFrame(out); o.to_csv(os.path.join(REV, "hotspot_allland_sensitivity.csv"),index=False)
print(o.round(2).to_string(index=False))
