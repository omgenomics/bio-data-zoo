# Contributing

## Suggest and contribute new formats

If you're interested in datasets for a file format or edge case that is not currently supported:

1. Create an issue on GitHub, describing the file format and scenario(s) to support.
2. Mention if you are interested in contributing the sample data.
3. If you want to contribute the data, wait for confirmation before you start working on it, and follow the instructions in the next section to submit a PR.


## Adding a new file

The general approach is: for each file format, we fetch sample files from public cloud buckets (logic in `Makefile`'s `init` target). Those are generally called `basic_[description].[format]`. Then, we modify those original files to simulate certain scenarios (logic in `generate_[format].sh`).

* If you're adding a scenario for an existing file format:
  * Modify `./src/generate_[format].sh` and add:

    ```bash
    # ------------------------------------------------------------------------------
    # Scenario
    # ------------------------------------------------------------------------------

    log "Creating [describe scenario]"
    # Add code here to generate that file. You can use $DIR_BASIC for inputs and $DIR_OUT for output
    head -n 10 "$DIR_BASIC" > "$DIR_OUT/new_scenario.fastq"
    # Write a command here within the "$()" that should return a non-empty string if the generation succeeded
    validate "$(diff "$DIR_BASIC" "$DIR_OUT/new_scenario.fastq")"
    ```

* If you're adding a new file format:
  * Create `./src/generate_[format].sh` and follow the structure of other similar files such as `generate_bed.sh`
  * Update the Makefile (we can work with you on that)
  * Update README to specify where the original files came from (we can work with you on that)

## Dev environment

Prerequisites to develop on this repo:

* `samtools`
* `bedtools`
* `bcftools`
* `tabix`
* `bgzip`
* `seqtk`
* `gshuf` (`brew install coreutils` on Mac)
* `zcat`
