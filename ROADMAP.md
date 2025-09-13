# 🚀 ChromR Project Roadmap

This roadmap outlines the planned improvements and future directions for **chromR**.  
Check the boxes as tasks are completed ✅

---

## 📌 Short-term goals (0–3 months)
- [ ] Improve documentation
  - [ ] Expand `README.md` with clear usage examples
  - [ ] Add installation instructions
  - [ ] Provide minimal working dataset for testing
- [ ] Improve outputs
  - [ ] Organize results into structured folders (`hmms/`, `search/`, `logs/`)
  - [ ] Standardize output formats (CSV, TSV, JSON)
- [ ] Add more HMMER functionality
  - [ ] Include `hmmscan` support
  - [ ] Add wrapper for `hmmpress`
- [ ] Testing & CI
  - [ ] Add unit tests (e.g., with `pytest` or `testthat`)
  - [ ] Integrate GitHub Actions for automated testing


---

## 📌 Mid-term goals (3–6 months)
- [ ] Add HPC/Cluster support
  - [ ] Generate SLURM/PBS submission scripts automatically
  - [ ] Provide `--slurm` flag for easy batch jobs
- [ ] Improve configurability
  - [ ] Support config files (YAML/JSON) to avoid long command lines
  - [ ] Add logging verbosity levels (`--quiet`, `--debug`)
- [ ] Performance improvements
  - [ ] Benchmark against large datasets
  - [ ] Parallelize searches where possible

---

## 📌 Long-term goals (6–12 months)
- [ ] Visualization
  - [ ] Plot HMM hit distributions across chromosomes
  - [ ] Generate summary heatmaps for hits
- [ ] Extend compatibility
  - [ ] Docker container for easy deployment
  - [ ] Package for Bioconda
- [ ] Community & contributions
  - [ ] Add contributing guidelines
  - [ ] Define coding style and standards
  - [ ] Encourage external contributions via issues & pull requests

---

## 📌 Stretch goals (future ideas)
- [ ] Web-based interface for running searches and visualizing results
- [ ] Integration with genome browsers (e.g., JBrowse, IGV)
- [ ] Machine learning approaches for classifying hits

---

✨ This roadmap is open for feedback. If you have ideas, feel free to open an issue or submit a PR!
