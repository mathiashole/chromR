# :honeybee: chromR

[![R](https://img.shields.io/badge/R-276DC3?style=for-the-badge&logo=R&logoColor=white&labelColor=101010)](https://www.r-project.org/about.html)
![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/mathiashole/chromR?style=for-the-badge&labelColor=101010&color=white)

`chromR` is an R script for customizing and visualizing gene or domain in chromosome or contig.

## :book: Features...

-  :bar_chart: Custom Visualization: Display genes, domains, or specific regions on chromosomes or contigs.
-  🎨 Color Customization: Define colors for categories using either hex color codes (e.g., #1f77b4) or standard R color names (e.g., black, orange).
-  :wrench: Keyword Filtering: Highlight specific genes, regions, or domains of interest using keywords.
-  :page_facing_up: GFF Support: Parse and plot data from GFF (General Feature Format) files.
-  :rocket: Command-Line Friendly: Run directly from the terminal or integrate into pipelines.
-  :computer: RStudio Compatible: Easy to modify, debug, or execute within RStudio.

## :wrench: Usage

To use the gene on chromosome plotting script, follow these steps:

1. Download or clone the script from the repository.
2. Ensure the required R packages are installed: `dplyr`, `readr` and `ggplot2`
3. **Option 1: Using the terminal**  
   Open a terminal, navigate to the folder containing the script, and run the script with the appropriate arguments.
4. **Option 2: Using RStudio**  
   Open the script directly in RStudio and run it from there, passing the arguments manually or setting them as variables within the script.

## :hammer: in progress ...

## :bulb: Bash Quick Examples

```{bash, eval = FALSE}
Rscript chromR.R -g </path/of/file.gff> -k <keyword> <keyword_asociated>
```

```{bash, eval = FALSE}
Rscript chromR.R -g </path/of/file.gff> -k <keyword1> <keyword1_asociated> <keyword2> <keyword2_asociated>
```

```{bash, eval = FALSE}
Rscript chromR.R -g </path/of/file.gff> -k <keyword> <keyword_asociated> -l </path/of/file_id>
```

## :sparkling_heart: Contributing

- :octocat: [Pull requests](https://github.com/mathiashole/chromR/pulls) and :star2: stars are always welcome.
- For major changes, please open an [issue](https://github.com/mathiashole/chromR/issues) first to discuss what you would like to change.
- Please make sure to update tests as appropriate.

## :mega: Contact

:bird: [Mathias](https://twitter.com/joaquinmangino)

## License
MIT &copy; [Mathias Mangino](https://github.com/mathiashole)