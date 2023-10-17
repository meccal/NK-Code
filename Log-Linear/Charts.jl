"""
LUCA MECCA
lmecca@london.edu
July 2023
"""

#################################################################################
################################### MP SHOCKS ###################################
#################################################################################
p = plot(layout = (5, 2), legend = false, size=(800,800))

plot!(p[1],IRF_output_gap_MP,
     ylim = (-0.4, 0),   # Set the y-axis boundaries
     xlim = (0, 15),   # Set the x-axis boundaries
     xticks = [0, 5, 10, 15],   # Set the desired x-axis ticks
     marker = (:circle, 4, :white),   # Use small empty circles as markers
     linecolor = :grey,   # Set the line color
     linewidth = 2,   # Set the line width
     legend = false,   # Hide the legend
     title = "Output gap",   # Set the title
     titlefontsize = 10,   # Adjust the size of the title
     grid = false   # Turn off the grid
)

plot!(p[2],IRF_inflation_MP,
     ylim = (-0.4, 0),   # Set the y-axis boundaries
     xlim = (0, 15),   # Set the x-axis boundaries
     xticks = [0, 5, 10, 15],   # Set the desired x-axis ticks
     marker = (:circle, 4, :white),   # Use small empty circles as markers
     linecolor = :grey,   # Set the line color
     linewidth = 2,   # Set the line width
     legend = false,   # Hide the legend
     title = "Inflation",   # Set the title
     titlefontsize = 10,   # Adjust the size of the title
     grid = false   # Turn off the grid
)

plot!(p[3],IRF_output_gap_MP,
     ylim = (-0.4, 0),   # Set the y-axis boundaries
     xlim = (0, 15),   # Set the x-axis boundaries
     xticks = [0, 5, 10, 15],   # Set the desired x-axis ticks
     marker = (:circle, 4, :white),   # Use small empty circles as markers
     linecolor = :grey,   # Set the line color
     linewidth = 2,   # Set the line width
     legend = false,   # Hide the legend
     title = "Output",   # Set the title
     titlefontsize = 10,   # Adjust the size of the title
     grid = false   # Turn off the grid
)

plot!(p[4],IRF_employment_MP,
     ylim = (-0.4, 0),   # Set the y-axis boundaries
     xlim = (0, 15),   # Set the x-axis boundaries
     xticks = [0, 5, 10, 15],   # Set the desired x-axis ticks
     marker = (:circle, 4, :white),   # Use small empty circles as markers
     linecolor = :grey,   # Set the line color
     linewidth = 2,   # Set the line width
     legend = false,   # Hide the legend
     title = "Employment",   # Set the title
     titlefontsize = 10,   # Adjust the size of the title
     grid = false   # Turn off the grid
)

plot!(p[5],IRF_real_wage_MP,
     ylim = (-2, 0),   # Set the y-axis boundaries
     xlim = (0, 15),   # Set the x-axis boundaries
     xticks = [0, 5, 10, 15],   # Set the desired x-axis ticks
     marker = (:circle, 4, :white),   # Use small empty circles as markers
     linecolor = :grey,   # Set the line color
     linewidth = 2,   # Set the line width
     legend = false,   # Hide the legend
     title = "Real wage",   # Set the title
     titlefontsize = 10,   # Adjust the size of the title
     grid = false   # Turn off the grid
)

plot!(p[6],IRF_price_level_MP,
     ylim = (-0.2, 0),   # Set the y-axis boundaries
     xlim = (0, 15),   # Set the x-axis boundaries
     xticks = [0, 5, 10, 15],   # Set the desired x-axis ticks
     marker = (:circle, 4, :white),   # Use small empty circles as markers
     linecolor = :grey,   # Set the line color
     linewidth = 2,   # Set the line width
     legend = false,   # Hide the legend
     title = "Price level",   # Set the title
     titlefontsize = 10,   # Adjust the size of the title
     grid = false   # Turn off the grid
)

plot!(p[7],IRF_nominal_rate_MP,
     ylim = (0, 0.4),   # Set the y-axis boundaries
     xlim = (0, 15),   # Set the x-axis boundaries
     xticks = [0, 5, 10, 15],   # Set the desired x-axis ticks
     marker = (:circle, 4, :white),   # Use small empty circles as markers
     linecolor = :grey,   # Set the line color
     linewidth = 2,   # Set the line width
     legend = false,   # Hide the legend
     title = "Nominal rate",   # Set the title
     titlefontsize = 10,   # Adjust the size of the title
     grid = false   # Turn off the grid
)

plot!(p[8],IRF_real_rate_MP,
     ylim = (0, 0.8),   # Set the y-axis boundaries
     xlim = (0, 15),   # Set the x-axis boundaries
     xticks = [0, 5, 10, 15],   # Set the desired x-axis ticks
     marker = (:circle, 4, :white),   # Use small empty circles as markers
     linecolor = :grey,   # Set the line color
     linewidth = 2,   # Set the line width
     legend = false,   # Hide the legend
     title = "Real rate",   # Set the title
     titlefontsize = 10,   # Adjust the size of the title
     grid = false   # Turn off the grid
)

plot!(p[9],IRF_money_supply_MP,
     ylim = (-0.8, 0),   # Set the y-axis boundaries
     xlim = (0, 15),   # Set the x-axis boundaries
     xticks = [0, 5, 10, 15],   # Set the desired x-axis ticks
     marker = (:circle, 4, :white),   # Use small empty circles as markers
     linecolor = :grey,   # Set the line color
     linewidth = 2,   # Set the line width
     legend = false,   # Hide the legend
     title = "Money supply",   # Set the title
     titlefontsize = 10,   # Adjust the size of the title
     grid = false   # Turn off the grid
)

plot!(p[10],MP_shocks,
     ylim = (0, 0.4),   # Set the y-axis boundaries
     xlim = (0, 15),   # Set the x-axis boundaries
     xticks = [0, 5, 10, 15],   # Set the desired x-axis ticks
     marker = (:circle, 4, :white),   # Use small empty circles as markers
     linecolor = :grey,   # Set the line color
     linewidth = 2,   # Set the line width
     legend = false,   # Hide the legend
     title = "v",   # Set the title
     titlefontsize = 10,   # Adjust the size of the title
     grid = false   # Turn off the grid
)

savefig("C:/Users/lmecca/OneDrive - London Business School/Research/Replications/NK/Output/Log-Linear/MP_shock.png")