using DelimitedFiles
using LinearAlgebra # I(m)

#max c' x
# s.t. Ax <= b, x >= 0

const TOL = 1e-10

struct ResultadoSimplex
  status::String
  value::Float64
  x::Vector{Float64}
end

ResultadoSimplex(status::String) = ResultadoSimplex(status, 0.0, Float64[])

function get_optimal_solution(tableau, n, m)
  optimal = tableau[end, end]
  
  x = zeros(n)

  for j in 1:n
    column = tableau[1:m, j]

    # so uma coluna com nao zeros
    if sum(abs.(column) .> TOL) == 1 
        idx = findfirst(x -> abs(x - 1.0) < TOL, column)
        if idx !== nothing
            x[j] = tableau[idx, end]  
        end
    end
  end

  return x, optimal
end

function read_problem(file::String)
  input_data = readdlm(file)

  m = Int(input_data[1, 1])
  n = Int(input_data[1, 2])

  c = input_data[2, 1:n] # da linha 2, na coluna 1 ate n
  A = input_data[3:end, 1:n]
  b = input_data[3:end, n+1]

  return n, m, A, b, c
end

function output_solution(file::String, result)
  open(file, "w") do f
    if result.status == "optimal"
      println(f, "otima")
      println(f, "valor otimo: ", result.value)
      println(f, "Solução:")
      for i in 1:length(result.x)
          println(f, "  x", i, " = ", result.x[i])
      end
    elseif result.status == "unlimited"
        println(f, "ilimitada")
        
    elseif result.status == "unviable"
        println(f, "inviável")
    end
  end
end

function generate_tableau(c, A, b)
  m, n = size(A)

  tableau = zeros(m+1, n + m + 1)

  # [A | I | b]
  # [c'| 0 | 0]

  tableau[1:m, 1:n] = A
  tableau[1:m, n+1:n+m] = I(m)
  tableau[1:m, end] = b

  tableau[end, 1:n] = -c

  return tableau
end

function choose_pivot_row(tableau, col)
  m = size(tableau, 1) - 1 # tamanho do A
  b = tableau[1:m, end] 
  column = tableau[1:m, col]

  ratios = fill(Inf, m)

  for i in 1:m
    if column[i] > TOL # tolerancia. pra n comparar com == 0
      ratios[i] = b[i] / column[i]
    end
  end

  min_ratio, row = findmin(ratios)

  return min_ratio < Inf ? row : -1 # se retornou -1
                                    # eh pq todos os anteriores
                                    # a[i, j] <= 0
                                    # logo, temos
end

function choose_pivot_col(tableau)
  # Regra: escolhe como pivot o menor valor, de menor indice

  last_row = tableau[end, 1:end-1] # aonde esta o C'

  min_val, col = findmin(last_row)

  return min_val < 0 ? col : -1 # se retornou -1
                                  # eh porque todos sao > 0
                                  # logo, eh otimo
end

function pivot!(tableau, row_pivot, col_pivot)
  m, n = size(tableau)

  pivot = tableau[row_pivot, col_pivot]
  # fazer pivot virar 1
  tableau[row_pivot, :] /= pivot

  # zerar o resto das linhas naquela coluna
  for i in 1:m
    if i != row_pivot
      factor = tableau[i, col_pivot]
      tableau[i, :] -= factor * tableau[row_pivot, :]
    end
  end
end

function simplex!(tableau, n, m)
  while true
    # var q entra
    pivot_col = choose_pivot_col(tableau)

    if pivot_col == -1 # Solucao otima 
      return get_optimal_solution(tableau, n, m)
    end

    # var que sai 
    pivot_row = choose_pivot_row(tableau, pivot_col)

    if pivot_row == -1
      error("unlimited")
    end

    pivot!(tableau, pivot_row, pivot_col)
  end
end

function main()
  if length(ARGS) < 2
    println("uso: julia simplex.jl <entrada> <saida>")
    exit(1)
  end

  input = ARGS[1]
  output = ARGS[2]

  n, m, A, b, c = read_problem(input)

  try
    tableau = generate_tableau(c, A, b)

    x, optimal_value = simplex!(tableau, n, m)

    result = ResultadoSimplex("optimal", optimal_value, x)

    output_solution(output, result)
  catch e
    println(string(e))
    if occursin("unlimited", string(e))
      output_solution(output, ResultadoSimplex("unlimited"))
    else
      output_solution(output, ResultadoSimplex("error"))
    end
  end
end

main()
