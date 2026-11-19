"""Re-download pscomppars from the NASA Exoplanet Archive TAP service (same query as the 2026-11-18 snapshot).
A new download will NOT reproduce md5 f2f70b24e67a9054940267e8081c81e8 (the archive is updated continuously);
the certified results use data/pscomppars_2026-11-18.csv. Usage: python fetch_catalog.py out.csv"""
import hashlib, sys, urllib.parse, urllib.request
Q = "select pl_name,hostname,discoverymethod,tran_flag,rv_flag,pl_controv_flag,pl_orbper,pl_orbpererr1,pl_orbpererr2,pl_orbsmax,pl_orbsmaxerr1,pl_orbsmaxerr2,pl_rade,pl_radeerr1,pl_radeerr2,pl_bmasse,pl_bmasseerr1,pl_bmasseerr2,pl_bmassprov,pl_orbeccen,pl_orbeccenerr1,pl_orbeccenerr2,pl_insol,pl_eqt,st_teff,st_tefferr1,st_tefferr2,st_rad,st_raderr1,st_raderr2,st_mass,st_masserr1,st_masserr2,st_lum,st_lumerr1,st_lumerr2,st_age,st_ageerr1,st_ageerr2,st_met,sy_dist,sy_jmag,sy_kmag,pl_rade_reflink,pl_bmasse_reflink from pscomppars"
url = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?" + urllib.parse.urlencode({"query": Q, "format": "csv"})
data = urllib.request.urlopen(url, timeout=180).read()
open(sys.argv[1] if len(sys.argv) > 1 else "pscomppars_new.csv", "wb").write(data)
print("md5", hashlib.md5(data).hexdigest())
