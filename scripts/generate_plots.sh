#!/bin/bash

# Wrapper around generate_plots.sage

if [[ -v TMUX ]]
then
    # Runs each given target in a separate process in a separate tmux pane
    for target in $@
    do
        tmux splitw -v -l 2 -t 0 sage generate_plots.sage $target
        echo "Launched $target"
    done
else
    # Fallback, not parallelized since output will be messed up
    for target in $@
    do
        sage generate_plots.sage $target
    done
fi
