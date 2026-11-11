# The explorer's help, written from its manual. Each page of docs/src/explorer
# becomes a topic of app/web/ui/js/help-content.js, and each heading a section
# under the anchor that Documenter gives it, so that a link between two pages
# of the manual works the same on the site and in the app. Links that leave
# the manual go to the site. The screenshots stay there: the app is in view.
#
#     julia docs/explorer_help.jl        # write help-content.js
#
# docs/make.jl takes the pages of the manual from here and stops when the file
# is not what they give.
module ExplorerHelp

using Markdown

const SITE = "https://ianlap.github.io/SigmaTau.jl/dev/"
const SRC = joinpath(@__DIR__, "src")
const OUT = normpath(@__DIR__, "..", "app", "web", "ui", "js", "help-content.js")

# The manual in order: the title of each page in the navigation, and the page.
const PAGES = [
    "Overview" => "explorer/index.md",
    "Install and run" => "explorer/install.md",
    "The side panel" => "explorer/panel.md",
    "Data" => "explorer/data.md",
    "Synthesis" => "explorer/synthesis.md",
    "Stability" => "explorer/stability.md",
    "Spectrum" => "explorer/spectrum.md",
    "Dynamic" => "explorer/dynamic.md",
    "Holdover" => "explorer/holdover.md",
    "Steering" => "explorer/steering.md",
    "Ensemble" => "explorer/ensemble.md",
    "Live" => "explorer/live.md",
    "Exports, sessions and reports" => "explorer/export.md",
    "Estimators and noise types" => "explorer/estimators.md",
]

# What the links of a page are resolved against.
struct Scope
    page::String                            # the page, under docs/src
    manual::Dict{String,Vector{String}}     # each page of the manual => the anchors of its headings
    docstrings::Dict{String,String}         # each documented name => the page of its docstring
    keys::Set{String}                       # the keys of the bibliography
end

parsepage(page) = Markdown.parse(read(joinpath(SRC, split(page, '/')...), String))
escape(s) = replace(s, '&' => "&amp;", '<' => "&lt;", '>' => "&gt;", '"' => "&quot;")

# ── Anchors and addresses ──

# The anchor of a heading, by the rule of Documenter.slugify.
function slug(s)
    s = replace(s, r"\s+" => "-")
    s = replace(s, "&" => "-and-")
    s = replace(s, r"[^\p{L}\p{P}\d\-]+" => "")
    return String(strip(replace(s, r"\-\-+" => "-"), '-'))
end

# The headings of a page with their anchors: the slug of the text, or of the
# name in `[text](@id name)`, and "-2", "-3" for a slug that the page repeats.
function headings(md)
    seen = Dict{String,Int}()
    return map(filter(b -> b isa Markdown.Header, md.content)) do h
        name = length(h.text) == 1 && h.text[1] isa Markdown.Link ? match(r"^@id (.+)$", h.text[1].url) : nothing
        text = name === nothing ? h.text : h.text[1].text
        id = slug(name === nothing ? sprint(Markdown.plain, Markdown.Paragraph(text)) : name[1])
        n = seen[id] = get(seen, id, 0) + 1
        (; level = typeof(h).parameters[1], text, anchor = n == 1 ? id : "$id-$n")
    end
end

# The page that a link on `from` names, as a path under docs/src.
function target(from, to)
    parts = split(from, '/')[1:end-1]
    for part in split(to, '/')
        part == ".." ? pop!(parts) : push!(parts, part)
    end
    return join(parts, '/')
end

topic(page) = splitext(basename(page))[1]
sitelink(page, anchor = "") =
    SITE * replace(page, r"(^|/)index\.md$" => s"\1", r"\.md$" => "/") * (isempty(anchor) ? "" : "#" * anchor)

# Where a link of the manual leads in the help: "#help/topic/anchor" inside
# the manual, a page of the site outside it.
function address(link, at)
    url = link.url
    startswith(url, "https://") && return url
    fail(why) = error("explorer help: $(at.page): the link ($url) $why")
    if url == "@ref"
        length(link.text) == 1 && link.text[1] isa Markdown.Code ||
            fail("is not to a docstring; write a link to a heading as page.md#Anchor")
        name = link.text[1].code
        haskey(at.docstrings, name) || fail("names nothing that an `@docs` block lists")
        return sitelink(at.docstrings[name], "SigmaTau." * name)
    elseif startswith(url, "@cite ")
        url[7:end] in at.keys || fail("names no record of refs.bib")
        return sitelink("bibliography.md", url[7:end])
    end
    if startswith(url, "@ref ")
        anchor = slug(strip(url[6:end], '"'))
        homes = [p for (p, anchors) in at.manual if anchor in anchors]
        length(homes) == 1 || fail("names no heading of the manual")
        page = homes[1]
    else
        path, anchor = occursin('#', url) ? split(url, '#'; limit = 2) : (url, "")
        page = isempty(path) ? at.page : target(at.page, path)
    end
    if haskey(at.manual, page)
        isempty(anchor) || anchor in at.manual[page] || fail("names no heading of $page")
        top = isempty(anchor) || anchor == first(at.manual[page])
        return "#help/" * topic(page) * (top ? "" : "/" * anchor)
    end
    endswith(page, ".md") && isfile(joinpath(SRC, split(page, '/')...)) || fail("names no page of the documentation")
    isempty(anchor) || any(h -> h.anchor == anchor, headings(parsepage(page))) || fail("names no heading of $page")
    return sitelink(page, anchor)
end

# ── Formulas ──

const THIN = " "

# A formula as HTML. It covers the TeX that the manual uses: letters, digits
# and signs, sub- and superscripts, `\sqrt`, `\frac`, `\sum`, `\ln` and `\,`.
# Anything else is an error, and a case to add here.
function tex(formula)
    cs = collect(formula)
    html, i = texrun(cs, 1, false, formula)
    i > length(cs) || error("explorer help: unmatched } in ``$formula``")
    return html
end

# From `i` to an unmatched "}" or the end: the HTML, and where it stopped.
# A sub- or superscript (`script`) takes no spaces around its operators.
function texrun(cs, i, script, formula)
    io = IOBuffer()
    operand = false     # whether what was written last can stand left of a binary operator
    while i <= length(cs) && cs[i] != '}'
        c = cs[i]
        i += 1
        if isspace(c)
        elseif c == '{'
            inner, i = texarg(cs, i - 1, script, formula)
            print(io, inner)
            operand = true
        elseif c == '_' || c == '^'
            tag = c == '_' ? "sub" : "sup"
            inner, i = texarg(cs, i, true, formula)
            print(io, "<$tag>", inner, "</$tag>")
            operand = true
        elseif c == '\\'
            j = something(findnext(!isletter, cs, i), length(cs) + 1)
            name = String(cs[i:max(i, j - 1)])      # a word, or one sign as in "\,"
            i += length(name)
            if name == ","
                print(io, THIN)
            elseif name == "sqrt"
                inner, i = texarg(cs, i, script, formula)
                print(io, "√<span class=\"rad\">", inner, "</span>")
                operand = true
            elseif name == "frac"
                above, i = texarg(cs, i, script, formula)
                below, i = texarg(cs, i, script, formula)
                print(io, "<sup>", above, "</sup>⁄<sub>", below, "</sub>")
                operand = true
            elseif name == "sum"
                print(io, "Σ")
                if get(cs, i, ' ') == '_'
                    inner, i = texarg(cs, i + 1, true, formula)
                    print(io, "<sub>", inner, "</sub>")
                end
                print(io, THIN)
                operand = false
            elseif name == "ln"
                next = get(cs, something(findnext(!isspace, cs, i), 0), ' ')
                print(io, operand ? THIN : "", "ln", next == '(' ? "" : THIN)
                operand = false
            else
                error("explorer help: ``$formula`` uses \\$name, which docs/explorer_help.jl does not convert")
            end
        elseif c in "=<>≤≥≈∝"
            print(io, script ? escape(string(c)) : " " * escape(string(c)) * " ")
            operand = false
        elseif c == '+' || c == '-'
            sign = c == '-' ? "−" : "+"
            print(io, operand && !script ? " $sign " : sign)
            operand = false
        else
            print(io, isascii(c) && isletter(c) ? "<i>$c</i>" : escape(string(c)))
            operand = isletter(c) || isdigit(c) || c == ')'
        end
    end
    return String(take!(io)), i
end

# The argument at `i`, a group in braces or one character: its HTML, and what follows it.
function texarg(cs, i, script, formula)
    get(cs, i, '\\') == '\\' && error("explorer help: ``$formula``: write the argument in braces")
    cs[i] == '{' || return first(texrun(cs[i:i], 1, script, formula)), i + 1
    inner, j = texrun(cs, i + 1, script, formula)
    j <= length(cs) || error("explorer help: unmatched { in ``$formula``")
    return inner, j + 1
end

# ── Markdown to HTML ──

inline(xs::AbstractVector, at) = join(inline(x, at) for x in xs)
inline(x::AbstractString, at) = escape(x)
inline(x::Markdown.Bold, at) = "<b>" * inline(x.text, at) * "</b>"
inline(x::Markdown.Italic, at) = "<em>" * inline(x.text, at) * "</em>"
inline(x::Markdown.Code, at) = "<code>" * escape(x.code) * "</code>"
inline(x::Markdown.LaTeX, at) = "<span class=\"math\">" * tex(x.formula) * "</span>"
inline(::Markdown.Image, at) = ""
inline(x, at) = error("explorer help: $(at.page) has a $(typeof(x)), which docs/explorer_help.jl does not convert")

function inline(x::Markdown.Link, at)
    body = inline(x.text, at)
    isempty(body) && return ""      # a link around a screenshot
    url = escape(address(x, at))
    return startswith(url, "#") ? "<a href=\"$url\">$body</a>" :
           "<a href=\"$url\" target=\"_blank\" rel=\"noopener\">$body</a>"
end

blocks(xs, at) = join(block(x, at) for x in xs)
block(x, at) = error("explorer help: $(at.page) has a $(typeof(x)), which docs/explorer_help.jl does not convert")

function block(x::Markdown.Paragraph, at)
    body = inline(x.content, at)
    return isempty(body) ? "" : "<p>$body</p>"
end

function block(x::Markdown.Code, at)
    startswith(x.language, "@") && error("explorer help: $(at.page) has a ```$(x.language) block")
    return "<pre><code>" * escape(x.code) * "</code></pre>"
end

function block(x::Markdown.List, at)
    tag = x.ordered == -1 ? "ul" : "ol"
    items = map(x.items) do item
        tight = length(item) == 1 && item[1] isa Markdown.Paragraph
        "<li>" * (tight ? inline(item[1].content, at) : blocks(item, at)) * "</li>"
    end
    return "<$tag>" * join(items) * "</$tag>"
end

function block(x::Markdown.Table, at)
    side = Dict(:l => "", :c => " class=\"c\"", :r => " class=\"r\"")
    row(cells, tag) = "<tr>" * join("<$tag$(side[a])>" * inline(c, at) * "</$tag>" for (c, a) in zip(cells, x.align)) * "</tr>"
    return "<table><thead>" * row(x.rows[1], "th") * "</thead><tbody>" *
           join(row(r, "td") for r in x.rows[2:end]) * "</tbody></table>"
end

block(x::Markdown.Admonition, at) =
    "<aside class=\"$(x.category)\"><h4>" * escape(x.title) * "</h4>" * blocks(x.content, at) * "</aside>"

# ── The file ──

# The topics of the help: for each page its sections, one a heading, each with
# its blocks of HTML. The first section is the page down to its first heading.
function topics()
    pages = last.(PAGES)
    manual = Dict(p => [h.anchor for h in headings(parsepage(p))] for p in pages)
    docstrings = Dict{String,String}()
    for (dir, _, files) in walkdir(SRC), file in filter(endswith(".md"), files)
        text = read(joinpath(dir, file), String)
        page = replace(relpath(joinpath(dir, file), SRC), '\\' => '/')
        for listing in eachmatch(r"^```@docs\s*$(.*?)^```"ms, text), name in split(listing[1])
            get!(docstrings, name, page)
        end
    end
    keys = Set(m[1] for m in eachmatch(r"^@\w+\{([^,\s]+),"m, read(joinpath(SRC, "refs.bib"), String)))
    return map(PAGES) do (nav, page)
        at = Scope(page, manual, docstrings, keys)
        md = parsepage(page)
        heads = headings(md)
        md.content[1] isa Markdown.Header{1} && count(h -> h.level == 1, heads) == 1 ||
            error("explorer help: $page must open with its one # heading")
        sections = []
        for b in md.content
            if b isa Markdown.Header
                h = heads[length(sections)+1]
                title = inline(h.text, at)
                push!(sections, (; id = isempty(sections) ? "" : h.anchor, h.level, title,
                                 html = ["<h$(h.level)>$title</h$(h.level)>"]))
            else
                html = block(b, at)
                isempty(html) || push!(sections[end].html, html)
            end
        end
        (; id = topic(page), nav, web = sitelink(page), sections)
    end
end

js(s) = '"' * replace(s, '\\' => "\\\\", '"' => "\\\"", '\n' => "\\n") * '"'

# The text of help-content.js.
function content()
    io = IOBuffer()
    println(io, "// The explorer's help: its manual, docs/src/explorer/*.md, as docs/explorer_help.jl")
    println(io, "// writes it. Change the pages and run that script; this file is not edited by hand.")
    println(io, "export const SITE = ", js(SITE), ";")
    println(io, "export const TOPICS = [")
    for t in topics()
        println(io, "  { id: ", js(t.id), ", nav: ", js(t.nav), ", web: ", js(t.web), ", sections: [")
        for s in t.sections
            println(io, "    { id: ", js(s.id), ", level: ", s.level, ", title: ", js(s.title), ", html: [")
            foreach(b -> println(io, "      ", js(b), ","), s.html)
            println(io, "    ] },")
        end
        println(io, "  ] },")
    end
    println(io, "];")
    return String(take!(io))
end

# Stop when help-content.js is not what the manual gives.
check() = isfile(OUT) && read(OUT, String) == content() ||
    error("app/web/ui/js/help-content.js is behind the manual: run `julia docs/explorer_help.jl` and commit the file")

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    write(ExplorerHelp.OUT, ExplorerHelp.content())
    println("wrote ", relpath(ExplorerHelp.OUT))
end
