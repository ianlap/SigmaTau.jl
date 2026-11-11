# Shell completion scripts. Static text (no analysis or file reads when the
# shell completes), generated from the option table: each command completes
# the options it accepts, and each option with fixed values completes them.
# FILE and the options that take a path keep the shell's own file completion.

cmd_options(cmd::String) = String[op.long for op in OPTS if cmd in op.cmds]
subcommands() = String[c for (c, _, _) in COMMAND_HELP]

# Values offered after an option; --devs takes a comma list of them.
completion_values(op::Opt) = op.long == "--devs" ? estimator_names() : op.choices
valued_options() = [op for op in OPTS if !isempty(completion_values(op))]
free_options() = String[op.long for op in OPTS if !isempty(op.metavar) && isempty(completion_values(op))]

spaced(v::Vector{String}) = join(v, " ")

function bash_script()
    io = IOBuffer()
    print(io, """
    # sigmatau bash completion.
    # Install:  sigmatau completion bash > /usr/local/etc/bash_completion.d/sigmatau
    #      or:  source <(sigmatau completion bash)
    _sigmatau() {
        local cur="\${COMP_WORDS[COMP_CWORD]}" prev="\${COMP_WORDS[COMP_CWORD-1]}" item
        case "\$prev" in
    """)
    for op in valued_options()
        if op.long == "--devs"
            print(io, "        --devs) item=\"\${cur##*,}\"; COMPREPLY=( \$(compgen -P \"\${cur%\"\$item\"}\" -W \"",
                  spaced(completion_values(op)), "\" -- \"\$item\") ); return ;;\n")
        else
            print(io, "        ", op.long, ") COMPREPLY=( \$(compgen -W \"", spaced(op.choices), "\" -- \"\$cur\") ); return ;;\n")
        end
    end
    print(io, "        ", join(free_options(), "|"), ") COMPREPLY=(); return ;;\n")
    print(io, """
        esac
        if [[ \$COMP_CWORD -eq 1 ]]; then
            COMPREPLY=( \$(compgen -W "$(spaced(subcommands()))" -- "\$cur") ); return
        fi
        local opts=""
        case "\${COMP_WORDS[1]}" in
            completion) [[ \$COMP_CWORD -eq 2 ]] && COMPREPLY=( \$(compgen -W "$(spaced(SHELLS))" -- "\$cur") ); return ;;
            help) [[ \$COMP_CWORD -eq 2 ]] && COMPREPLY=( \$(compgen -W "$(spaced(subcommands()))" -- "\$cur") ); return ;;
    """)
    for c in COMMANDS
        print(io, "        ", c, ") opts=\"", spaced(cmd_options(c)), "\" ;;\n")
    end
    print(io, """
        esac
        if [[ "\$cur" == -* ]]; then
            COMPREPLY=( \$(compgen -W "\$opts" -- "\$cur") )
        else
            COMPREPLY=()    # FILE: the shell's file completion (-o default)
        fi
    }
    complete -o bashdefault -o default -F _sigmatau sigmatau
    """)
    return String(take!(io))
end

function zsh_script()
    io = IOBuffer()
    print(io, """
    #compdef sigmatau
    # sigmatau zsh completion.
    # Install:  sigmatau completion zsh > "\${fpath[1]}/_sigmatau"   (then restart the shell)
    _sigmatau() {
      local curcontext="\$curcontext" state line
      typeset -A opt_args
      _arguments -C '1: :->cmd' '*:: :->args'
      case \$state in
        cmd) _values 'command' $(spaced(subcommands())) ;;
        args)
          case \$words[1] in
            completion) _values 'shell' $(spaced(SHELLS)) ;;
            help) _values 'command' $(spaced(subcommands())) ;;
            *)
              case \$words[CURRENT-1] in
    """)
    for op in valued_options()
        sep = op.long == "--devs" ? "-s , " : ""
        print(io, "            ", op.long, ") _values ", sep, "'", op.metavar, "' ", spaced(completion_values(op)), " ;;\n")
    end
    print(io, "            ", join(free_options(), "|"), ") _files ;;\n")
    print(io, """
                *)
                  if [[ \$PREFIX == -* ]]; then
                    local -a opts
                    case \$words[1] in
    """)
    for c in COMMANDS
        print(io, "                  ", c, ") opts=(", spaced(cmd_options(c)), ") ;;\n")
    end
    print(io, """
                    esac
                    compadd -a opts
                  else
                    _files
                  fi ;;
              esac ;;
          esac ;;
      esac
    }
    _sigmatau "\$@"
    """)
    return String(take!(io))
end

fish_quote(s::String) = "'" * replace(s, "\\" => "\\\\", "'" => "\\'") * "'"

function fish_script()
    io = IOBuffer()
    print(io, """
    # sigmatau fish completion.
    # Install:  sigmatau completion fish > ~/.config/fish/completions/sigmatau.fish
    function __sigmatau_devs
        printf '%s\\n' $(spaced(estimator_names()))
    end
    complete -c sigmatau -f -n '__fish_use_subcommand' -a '$(spaced(subcommands()))'
    complete -c sigmatau -f -n '__fish_seen_subcommand_from completion' -a '$(spaced(SHELLS))'
    complete -c sigmatau -f -n '__fish_seen_subcommand_from help' -a '$(spaced(subcommands()))'
    """)
    for op in OPTS
        print(io, "complete -c sigmatau -n '__fish_seen_subcommand_from ", spaced(op.cmds), "' -l ", op.long[3:end])
        if op.long == "--devs"
            print(io, " -x -a '(__fish_complete_list , __sigmatau_devs)'")
        elseif !isempty(op.choices)
            print(io, " -x -a '", spaced(op.choices), "'")
        elseif op.path
            print(io, " -r -F")
        elseif !isempty(op.metavar)
            print(io, " -x")
        end
        print(io, " -d ", fish_quote(op.help), "\n")
    end
    return String(take!(io))
end

ps_table(pairs) = join(("        '" * k * "' = '" * v * "'" for (k, v) in pairs), "\n")

function powershell_script()
    values = [(op.long, spaced(completion_values(op))) for op in valued_options()]
    opts = [(c, spaced(cmd_options(c))) for c in COMMANDS]
    free = join(("'" * f * "'" for f in free_options()), ", ")
    return """
    # sigmatau PowerShell completion.
    # Install (add to \$PROFILE):  sigmatau completion powershell | Out-String | Invoke-Expression
    Register-ArgumentCompleter -Native -CommandName sigmatau -ScriptBlock {
        param(\$wordToComplete, \$commandAst, \$cursorPosition)
        \$values = @{
    $(ps_table(values))
        }
        \$options = @{
    $(ps_table(opts))
        }
        \$elements = @(\$commandAst.CommandElements | ForEach-Object { \$_.ToString() })
        \$sub  = if (\$elements.Count -ge 2) { \$elements[1] } else { '' }
        \$prev = if (\$wordToComplete) { if (\$elements.Count -ge 2) { \$elements[-2] } else { '' } } else { \$elements[-1] }
        # A comma list (--devs) completes its last item.
        \$head = if (\$prev -eq '--devs') { \$wordToComplete -replace '[^,]*\$', '' } else { '' }
        function script:_res(\$items) {
            \$items -split ' ' | ForEach-Object { \$head + \$_ } | Where-Object { \$_ -like "\$wordToComplete*" } | ForEach-Object {
                [System.Management.Automation.CompletionResult]::new(\$_, \$_, 'ParameterValue', \$_)
            }
        }
        if (\$values.ContainsKey(\$prev)) { return _res \$values[\$prev] }
        if (\$prev -in @($free)) { return }
        if (\$elements.Count -le 2 -and -not \$wordToComplete.StartsWith('-')) {
            if (\$elements.Count -le 1 -or \$elements[1] -eq \$wordToComplete) { return _res '$(spaced(subcommands()))' }
        }
        if (\$sub -eq 'completion') { return _res '$(spaced(SHELLS))' }
        if (\$sub -eq 'help') { return _res '$(spaced(subcommands()))' }
        if (\$wordToComplete.StartsWith('-') -and \$options.ContainsKey(\$sub)) { return _res \$options[\$sub] }
    }
    """
end

function completion_script(sh::String)
    sh == "bash" && return bash_script()
    sh == "zsh" && return zsh_script()
    sh == "fish" && return fish_script()
    (sh == "powershell" || sh == "pwsh") && return powershell_script()
    return ""
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
