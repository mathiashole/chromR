# :honeybee: chromR

[![R](https://img.shields.io/badge/R-276DC3?style=for-the-badge&logo=R&logoColor=white&labelColor=101010)](https://www.r-project.org/about.html)
![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/mathiashole/chromR?style=for-the-badge&labelColor=101010&color=white)

 :honeybee: `chromR` is an R script designed for visualizing gene or domain positions on chromosomes or contigs. It allows users to customize their plots using keywords, file-based inputs, and color options, making it ideal for genomic data exploration.

## :book: Features

-  📊 Custom Visualization: Display genes, domains, or specific regions on chromosomes or contigs.
-  🎨 Color Customization: Define colors for categories using either hex color codes (e.g., #1f77b4) or standard R color names (e.g., black, orange).
-  🛠️ Keyword Filtering: Highlight specific genes, regions, or domains of interest using keywords.
-  📄 GFF Support: Parse and plot data from GFF (General Feature Format) files.
-  🚀 Command-Line Friendly: Run directly from the terminal or integrate into pipelines.
-  🖥️ RStudio Compatible: Easy to modify, debug, or execute within RStudio.

## :wrench: Usage

### Requirements:

Ensure the following R packages are installed: `dplyr`, `readr`, `ggplot2`

### Execution Options:

#### 1. Using the Terminal
-  Navigate to the folder containing `chromR.R` and execute the script with the necessary arguments:

```{bash, eval = FALSE}
Rscript chromR.R --gff_file </path/to/file.gff> --keywords <keyword1> <keyword2> --colors <color1>
```

#### 2. Using RStudio

-  Open chromR.R in RStudio.
-  Set the arguments manually in the script or pass them interactively.

## :gear: Arguments

| Argument            | Description                               | Example                     |
|---------------------|-------------------------------------------|-----------------------------|
| `--gff_file` / `-g` | Path to the GFF file.                    | `--gff_file data.gff`       |
| `--keywords` / `-k` | Keywords to highlight specific features. | `--keywords gene1 CDS`      |
| `--colors` / `-c`   | Colors for keywords (names or hex).      | `--colors black orange`     |
| `--layout` / `-l`   | Optional layout input (e.g., ID file).   | `--layout data_ids.txt`     |
| `--strict` / `-s`   | Optional shows only the chromosomes that contain the keywords | `--strict`               |
| `--line_plot` / `-lp`| Optional line plot gene.                | `--line_plot`               |

:bulb: Note: Colors must match the number of keywords provided. They can be standard R color names or hex codes. In the case of hexadecimal colors, they must be put in quotes (e.g. "#F4A582").


## :hammer: in progress ...

## :bulb: Bash Quick Examples

Example 1: Basic Plot

Highlight gene1 and CDS in orange:

```{bash, eval = FALSE}
Rscript chromR.R --gff_file data.gff --keywords gene1 CDS --colors orange
```

```{bash, eval = FALSE}
Rscript chromR.R -g </path/of/file.gff> -k <keyword1> <keyword1_asociated> <keyword2> <keyword2_asociated>
```

```{bash, eval = FALSE}
Rscript chromR.R -g </path/of/file.gff> -k <keyword> <keyword_asociated> -l </path/of/file_id>
```

## :sparkling_heart: Contributing

We welcome contributions!

- :octocat: [Pull requests](https://github.com/mathiashole/chromR/pulls) and :star2: stars are always welcome.
- For major changes, please open an [issue](https://github.com/mathiashole/chromR/issues) first to discuss what you would like to change.
- Please make sure to update tests as appropriate.

## :mega: Contact

:bird: [Mathias](https://twitter.com/joaquinmangino)

## :scroll: License
MIT &copy; [Mathias Mangino](https://github.com/mathiashole)