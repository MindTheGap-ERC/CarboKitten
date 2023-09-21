using Documenter, CarboKitten

module Entangled
    using DataStructures: DefaultDict

    function transpile_md(src)
        counts = DefaultDict(0)
        Channel{String}() do ch
            for line in src
                if (m = match(r"``` *{[^#}]*#([a-zA-Z0-9\-_]+)[^}]*\}", line)) !== nothing
                    term = counts[m[1]] == 0 ? "≣" : "⊞"
                    put!(ch, "```@raw html")
                    put!(ch, "<div class=\"noweb-label\">⪡" * m[1] * "⪢" * term * "</div>")
                    put!(ch, "```")
                    put!(ch, line)
                    counts[m[1]] += 1
                elseif (m = match(r"``` *{[^}]*file=([a-zA-Z0-9\-_\.\/\\]+)[^}]*}", line)) !== nothing
                    put!(ch, "```@raw html")
                    put!(ch, "<div class=\"noweb-label\">file:<i>" * m[1] * "</i></div>")
                    put!(ch, "```")
                    put!(ch, line)
                else
                    put!(ch, line)
                end
            end
        end
    end

    function transpile_file(src, target_path)
        mkpath(joinpath(target_path, dirname(src)))
        content = open(readlines, src, "r")
        open(joinpath(target_path, basename(src)), "w") do fout
            join(fout, transpile_md(content), "\n")
        end
    end
end

is_markdown(path) = splitext(path)[2] == ".md"
sources = filter(is_markdown, readdir(joinpath(@__DIR__, "src"), join=true))
path = mktempdir()
Entangled.transpile_file.(sources, path)

makedocs(
    source=path,
    sitename="CarboKitten",
    pages = [
        "Bosscher and Schlager 1992" => "bosscher-1992.md",
        "CarboCAT" => [
            "summary" => "carbocat.md",
            "cellular automaton" => "carbocat-ca.md",
            "model with ca and production" => "ca-with-production.md",
            "sediment transport" => "carbocat-transport.md"
        ],
        "Algorithms" => [
            "stencils" => "stencils.md",
            "utility" => "utility.md"
        ]
    ])
