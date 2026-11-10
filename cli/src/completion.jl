# Shell completion scripts. Static text (no analysis or file reads when the
# shell completes), with every vocabulary filled in from the option table.
# FILE and --out keep the shell's own file completion.

cmd_options(cmd::String) = String[op.long for op in OPTS if cmd in op.cmds]
value_options() = String[op.long for op in OPTS if !isempty(op.metavar) && isempty(op.choices) && !op.path]
subcommands() = String[c for (c, _, _) in COMMAND_HELP]

function fill_template(tpl::String)
    s = tpl
    s = replace(s, "@@SUBCMDS@@" => join(subcommands(), " "))
    s = replace(s, "@@SHELLS@@" => join(SHELLS, " "))
    s = replace(s, "@@DEVS@@" => join(estimator_names(), " "))
    s = replace(s, "@@DYN@@" => join(DYNAMIC_DEVS, " "))
    s = replace(s, "@@GRIDS@@" => join(GRIDS, " "))
    s = replace(s, "@@FMT@@" => join(FORMATS, " "))
    s = replace(s, "@@DETREND@@" => join(DETRENDS, " "))
    s = replace(s, "@@DRIFT@@" => join(DRIFTS, " "))
    s = replace(s, "@@MODELS@@" => join(DRIFT_MODELS, " "))
    s = replace(s, "@@VALUEOPTS@@" => join(value_options(), "|"))
    for c in COMMANDS
        s = replace(s, "@@OPTS_" * c * "@@" => join(cmd_options(c), " "))
    end
    return s
end

const BASH = raw"""
# sigmatau bash completion.
# Install:  sigmatau completion bash > /usr/local/etc/bash_completion.d/sigmatau
#      or:  source <(sigmatau completion bash)
_sigmatau() {
    local cur prev
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"
    case "$prev" in
        --devs)    COMPREPLY=( $(compgen -W "@@DEVS@@" -- "${cur##*,}") ); return ;;
        --dev)     COMPREPLY=( $(compgen -W "@@DYN@@" -- "$cur") ); return ;;
        --grid)    COMPREPLY=( $(compgen -W "@@GRIDS@@" -- "$cur") ); return ;;
        --format)  COMPREPLY=( $(compgen -W "@@FMT@@" -- "$cur") ); return ;;
        --detrend) COMPREPLY=( $(compgen -W "@@DETREND@@" -- "$cur") ); return ;;
        --remove-drift) COMPREPLY=( $(compgen -W "@@DRIFT@@" -- "$cur") ); return ;;
        --model)   COMPREPLY=( $(compgen -W "@@MODELS@@" -- "$cur") ); return ;;
        --out)     COMPREPLY=(); return ;;
        @@VALUEOPTS@@) COMPREPLY=(); return ;;
    esac
    if [[ $COMP_CWORD -eq 1 ]]; then
        COMPREPLY=( $(compgen -W "@@SUBCMDS@@" -- "$cur") ); return
    fi
    local opts=""
    case "${COMP_WORDS[1]}" in
        completion) [[ $COMP_CWORD -eq 2 ]] && COMPREPLY=( $(compgen -W "@@SHELLS@@" -- "$cur") ); return ;;
        help)       [[ $COMP_CWORD -eq 2 ]] && COMPREPLY=( $(compgen -W "@@SUBCMDS@@" -- "$cur") ); return ;;
        analyze)    opts="@@OPTS_analyze@@" ;;
        check)      opts="@@OPTS_check@@" ;;
        spectrum)   opts="@@OPTS_spectrum@@" ;;
        dynamic)    opts="@@OPTS_dynamic@@" ;;
        drift)      opts="@@OPTS_drift@@" ;;
        holdover)   opts="@@OPTS_holdover@@" ;;
        report)     opts="@@OPTS_report@@" ;;
        live)       opts="@@OPTS_live@@" ;;
    esac
    if [[ "$cur" == -* ]]; then
        COMPREPLY=( $(compgen -W "$opts" -- "$cur") )
    else
        COMPREPLY=()    # FILE: the shell's file completion (-o default)
    fi
}
complete -o bashdefault -o default -F _sigmatau sigmatau
"""

const ZSH = raw"""
#compdef sigmatau
# sigmatau zsh completion.
# Install:  sigmatau completion zsh > "${fpath[1]}/_sigmatau"   (then restart the shell)
_sigmatau() {
  local curcontext="$curcontext" state line
  typeset -A opt_args
  _arguments -C '1: :->cmd' '*:: :->args'
  case $state in
    cmd)  _values 'command' @@SUBCMDS@@ ;;
    args)
      case $words[1] in
        completion) _values 'shell' @@SHELLS@@ ;;
        help)       _values 'command' @@SUBCMDS@@ ;;
        analyze|check|spectrum|dynamic|drift|holdover|report|live)
          case $words[CURRENT-1] in
            --devs)    _values -s , 'deviation' @@DEVS@@ ;;
            --dev)     _values 'deviation' @@DYN@@ ;;
            --grid)    _values 'grid' @@GRIDS@@ ;;
            --format)  _values 'format' @@FMT@@ ;;
            --detrend) _values 'mode' @@DETREND@@ ;;
            --remove-drift) _values 'mode' @@DRIFT@@ ;;
            --model)   _values 'model' @@MODELS@@ ;;
            --out)     _files ;;
            *)
              if [[ $PREFIX == -* ]]; then
                local -a opts
                case $words[1] in
                  analyze)  opts=(@@OPTS_analyze@@) ;;
                  check)    opts=(@@OPTS_check@@) ;;
                  spectrum) opts=(@@OPTS_spectrum@@) ;;
                  dynamic)  opts=(@@OPTS_dynamic@@) ;;
                  drift)    opts=(@@OPTS_drift@@) ;;
                  holdover) opts=(@@OPTS_holdover@@) ;;
                  report)   opts=(@@OPTS_report@@) ;;
                  live)     opts=(@@OPTS_live@@) ;;
                esac
                compadd -a opts
              else
                _files
              fi ;;
          esac ;;
      esac ;;
  esac
}
_sigmatau "$@"
"""

const FISH = raw"""
# sigmatau fish completion.
# Install:  sigmatau completion fish > ~/.config/fish/completions/sigmatau.fish
complete -c sigmatau -f -n '__fish_use_subcommand' -a '@@SUBCMDS@@'
complete -c sigmatau -f -n '__fish_seen_subcommand_from completion' -a '@@SHELLS@@'
complete -c sigmatau -f -n '__fish_seen_subcommand_from help' -a '@@SUBCMDS@@'
complete -c sigmatau -x -n '__fish_seen_subcommand_from analyze report live' -l devs -a '@@DEVS@@' -d 'comma list of deviations'
complete -c sigmatau -x -n '__fish_seen_subcommand_from dynamic' -l dev -a '@@DYN@@' -d 'dynamic deviation'
complete -c sigmatau -x -n '__fish_seen_subcommand_from analyze dynamic report live' -l grid -a '@@GRIDS@@' -d 'τ grid'
complete -c sigmatau -x -n '__fish_seen_subcommand_from drift' -l model -a '@@MODELS@@' -d 'drift model'
complete -c sigmatau -x -n '__fish_seen_subcommand_from analyze check spectrum dynamic drift holdover' -l format -a '@@FMT@@' -d 'output format'
complete -c sigmatau -x -n '__fish_seen_subcommand_from analyze check spectrum dynamic drift holdover report' -l detrend -a '@@DETREND@@' -d 'detrend mode'
complete -c sigmatau -x -n '__fish_seen_subcommand_from analyze check spectrum dynamic drift holdover report' -l remove-drift -a '@@DRIFT@@' -d 'drift to fit and remove'
complete -c sigmatau -r -F -n '__fish_seen_subcommand_from analyze check spectrum dynamic drift holdover report' -l out -d 'write to PATH'
"""

const POWERSHELL = raw"""
# sigmatau PowerShell completion.
# Install (add to $PROFILE):  sigmatau completion powershell | Out-String | Invoke-Expression
Register-ArgumentCompleter -Native -CommandName sigmatau -ScriptBlock {
    param($wordToComplete, $commandAst, $cursorPosition)
    $elements = @($commandAst.CommandElements | ForEach-Object { $_.ToString() })
    $sub  = if ($elements.Count -ge 2) { $elements[1] } else { '' }
    $prev = if ($wordToComplete) { if ($elements.Count -ge 2) { $elements[-2] } else { '' } } else { $elements[-1] }
    function script:_res($items) {
        $items | Where-Object { $_ -like "$wordToComplete*" } | ForEach-Object {
            [System.Management.Automation.CompletionResult]::new($_, $_, 'ParameterValue', $_)
        }
    }
    switch ($prev) {
        '--devs'    { return _res ('@@DEVS@@' -split ' ') }
        '--dev'     { return _res ('@@DYN@@' -split ' ') }
        '--grid'    { return _res ('@@GRIDS@@' -split ' ') }
        '--format'  { return _res ('@@FMT@@' -split ' ') }
        '--detrend' { return _res ('@@DETREND@@' -split ' ') }
        '--remove-drift' { return _res ('@@DRIFT@@' -split ' ') }
        '--model'   { return _res ('@@MODELS@@' -split ' ') }
        '--out'     { return }
    }
    if ($elements.Count -le 2 -and -not $wordToComplete.StartsWith('-')) {
        if ($elements.Count -le 1 -or $elements[1] -eq $wordToComplete) { return _res ('@@SUBCMDS@@' -split ' ') }
    }
    if ($sub -eq 'completion') { return _res ('@@SHELLS@@' -split ' ') }
    if ($sub -eq 'help') { return _res ('@@SUBCMDS@@' -split ' ') }
    if ($wordToComplete.StartsWith('-')) {
        switch ($sub) {
            'analyze'  { return _res ('@@OPTS_analyze@@' -split ' ') }
            'check'    { return _res ('@@OPTS_check@@' -split ' ') }
            'spectrum' { return _res ('@@OPTS_spectrum@@' -split ' ') }
            'dynamic'  { return _res ('@@OPTS_dynamic@@' -split ' ') }
            'drift'    { return _res ('@@OPTS_drift@@' -split ' ') }
            'holdover' { return _res ('@@OPTS_holdover@@' -split ' ') }
            'report'   { return _res ('@@OPTS_report@@' -split ' ') }
            'live'     { return _res ('@@OPTS_live@@' -split ' ') }
        }
    }
}
"""

function completion_script(sh::String)
    sh == "bash" && return fill_template(BASH)
    sh == "zsh" && return fill_template(ZSH)
    sh == "fish" && return fill_template(FISH) * fish_options()
    (sh == "powershell" || sh == "pwsh") && return fill_template(POWERSHELL)
    return ""
end

# fish: one line per option, scoped to the commands that take it.
function fish_options()
    io = IOBuffer()
    for op in OPTS
        # --help, --out and the options with fixed values are in FISH above.
        (op.long == "--help" || op.path || !isempty(op.choices)) && continue
        print(io, "complete -c sigmatau ")
        isempty(op.metavar) || print(io, op.path ? "-r -F " : "-x ")
        print(io, "-n '__fish_seen_subcommand_from ", join(op.cmds, " "), "' -l ", op.long[3:end])
        print(io, " -d '", replace(op.help, "'" => "\\'"), "'\n")
    end
    return String(take!(io))
end

function run_completion(args::Vector{String}, out::IO, err::IO)
    if length(args) < 2
        print(err, "sigmatau: completion needs a shell: ", join(SHELLS, ", "), "\n")
        return EXIT_USAGE
    end
    sh = lowercase(args[2])
    if sh == "-h" || sh == "--help"
        print(out, "Usage: sigmatau completion SHELL\n\nShells: ", join(SHELLS, ", "), "\n")
        return EXIT_OK
    end
    script = completion_script(sh)
    if isempty(script)
        print(err, "sigmatau: unknown shell \"", args[2], "\"; valid: ", join(SHELLS, ", "), "\n")
        return EXIT_USAGE
    end
    print(out, script)
    return EXIT_OK
end
