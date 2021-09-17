library(readr)

## ---- Read the file with the data you want to expand ----
cities <- read_csv("cities.csv")


## ---- Setup the features of the application ----

"Input the min/max coordinates you are working with here"
xmin <- -52
xmax <- -40.5
ymin <- -5.8
ymax <- 1

"Defines how many nodes you want to create inside your map. 
It can improve resolution, but hurt performance."
nodes <- 5

"Since we are completing the data of the map, this value is necessary
in order to complete the empty data on the map."
tx_padrao <- 30

"Just for some time measures, don't bother."
error = 1.5

"Choose 0 for default radial/square spread or 1 for diagonal only spread"
spread_option <- 0

"Choose the number of iterations. This represent how deep in you map
you want to spread your data over. Not bigger than your map width's 
or heigth's"
n_iterations <- 5


"The degree value setup how much each iteration upgrades spread data on
the next value adjacent."
degree <- 2


#Here the code starts
limits <- list(xmin, xmax, ymin, ymax)

# Generate a grid:
mtr = matrix(data = tx_padrao, nrow = round((ymax - ymin)*nodes), ncol = round((xmax-xmin)*nodes))

# Transform the coordinates into indexes for the matrix
transform_x <- function(x, limits, mtr){
  new_x <- round(((x - limits[[1]])/(limits[[2]] - limits[[1]]))*ncol(mtr), digits = 0)
  return(as.integer(new_x))
}

transform_y <- function(y, limits, mtr){
  new_y <- round(((y - limits[[3]])/(limits[[4]] - limits[[3]]))*nrow(mtr), digits = 0)
  return(as.integer(nrow(mtr) - new_y))
}

# Define the position in the matrix
pos <- function(x, y, mat){
  return((x-1)*nrow(mat) + y)
}

# Create a list of x, y coordinates for the cities
cities_coord <- list()
for (i in 1:nrow(cities)){
  new_x <- transform_x(cities$x[i], limits, mtr)
  new_y <- transform_y(cities$y[i], limits, mtr)
  mtr[pos(new_x, new_y, mtr)] <- round(cities$tc_ac[i], digits = 0)
  cities_coord <- append(cities_coord, list(new_x, new_y))
}

# Default modifiers -- Don't change them
linear_modifier <- list(1, 0, -1, 0, 0, 1, 0, -1)
linear_modifier2 <- Map('*', linear_modifier, 2)

##---- Spread Functions ----
linear_spread <- function(coords, matx, number_iterations, degree){
  
  linear_modifier_base <- list(1, 0, -1, 0, 0, 1, 0, -1) # D, E, C, B
  diagonal_modifier_base <- list(1, 1, -1, 1, -1, -1, 1, -1) # NE, NW, SW, SE
  
  for (i in 1:number_iterations){
    
    linear_modifier <- Map('*', linear_modifier_base, i)
    diagonal_modifier <- Map('*', diagonal_modifier_base, i)
    
    for (n in 1:(length(coords)/2)){
      
      coords_linear <- list()
      coords_diag <- list()
      
      new_x = coords[[n*2-1]]
      new_y = coords[[n*2]]
      position <- pos(new_x, new_y, matx)
      value <- matx[position]
      spread_line <- value - degree*i
      spread_diagonal <- value - degree*1.4*i
      for (cor in 1:4){
        coords_linear <- append(coords_linear, new_x)
        coords_linear <- append(coords_linear, new_y)
        coords_diag <- append(coords_diag, new_x)
        coords_diag <- append(coords_diag, new_y)
      }
      new_coords_linear <- Map( "+" , coords_linear, linear_modifier)
      new_coords_diagonal <- Map( "+" , coords_diag, diagonal_modifier)
      for (t in 1:4){
        linear_position <- pos(new_coords_linear[[t*2-1]], new_coords_linear[[t*2]], matx)
        if (matx[linear_position] == 30){
          matx[linear_position] <- spread_line
        }
        else {
          matx[linear_position] <- (spread_line + matx[linear_position])/2
        }
        diagonal_position <- pos(new_coords_diagonal[[t*2-1]], new_coords_diagonal[[t*2]], matx)
        if (matx[diagonal_position] == 30){
          matx[diagonal_position] <- spread_diagonal
        }
        else{
          matx[diagonal_position] <- (spread_diagonal + matx[diagonal_position])/2
        }
      }
    }
  }
  
  return(matx)
}

new_value <- function(x, y, old_value, degree){
  return(old_value - degree * sqrt(x ** 2 + y ** 2))
}

square_spread <- function(coordinates, matx, number_iterations, degree){
  
  for (i in 1:number_iterations){
    
    max_distance <- (2 * i + 1) / 2
    modifiers <- list()
    
    for (mod_x in 1 :( 2 * i + 1)){
      for (mod_y in 1 :( 2 * i + 1)){
        modifier_x <- mod_x - 2 * i
        modifier_y <- mod_y - 2 * i
        modifiers <- append(modifiers, modifier_x)
        modifiers <- append(modifiers, modifier_y)
      }
    }
    
    for (n in 1:(length(coordinates)/2)){

      x_coord = coordinates[[n*2-1]]
      y_coord = coordinates[[n*2]]
      position <- pos(x_coord, y_coord, matx)
      value <- matx[position]
      base_array = list()
      
      for (base_points in 1:(length(modifiers)/2)){
        base_array <- append(base_array, x_coord)
        base_array <- append(base_array, y_coord)
      }
      
      new_coords <- Map('+', base_array, modifiers)
      
      for (coord in 1:(length(new_coords)/2)){
        x <- new_coords[[coord * 2 - 1]]
        y <- new_coords[[coord * 2]]
        x_modifier <- modifiers[[coord * 2 - 1]]
        y_modifier <- modifiers[[coord * 2]]
        
        new_value <- new_value(x_modifier, y_modifier, value, degree)
        position <- pos(x, y, matx)
        
        if (matx[position] <= 30 && new_value > 30){
          matx[position] <- round(new_value, digits = 0)
        }
        else if (matx[position] > 30 && new_value > 30) {
          matx[position] <- round((matx[position] + new_value)/2, digits = 0)
        }
      }
    }
  }
  return(matx)
}


## ---- Now, the data really starts to be generated ----

if (spread_option == 0){
  new_mat <- square_spread(cities_coord, mtr, n_iterations, degree)  
} else if (spread_option == 0){
  new_mat <- square_spread(cities_coord, mtr, n_iterations, degree)
} else {
  error("Not possible.")
}

# Function to convert matrix indexes back to coordinates
trans_coord_x <- function(x, limits, mat){
  a = limits[[1]]
  b = limits[[2]]
  c = 1
  d = ncol(mat)
  X_map = (x*(b - a) + d*a - b*c)/(d-c)
}

trans_coord_y <- function(y, limits, mat){
  a = limits[[3]]
  b = limits[[4]]
  c = nrow(mat)
  d = 1
  Y_map = (y*(b - a) + d*a - b*c)/(d-c)
}

# Some metrics for time measure
linhas <- nrow(new_mat)
colunas <- ncol(new_mat)
total_lines <- linhas * colunas * mean(new_mat)

print_UI <- function(time_passed, expected_time){
  sprintf('Time passed: %.2f minutes ... Expectation: %.2f minutes', time_passed, expected_time)
}

# Final script to create the data
for (l in 1:linhas){
  if (l == 1){
    start_time <- Sys.time()
  }
  for (c in 1:colunas){
    start_loop_time <- Sys.time()
    x1 = trans_coord_x(c, limits, new_mat)
    y1 = trans_coord_y(l, limits, new_mat)
    
    pos_value <- pos(c, l, new_mat)
    n <- round(new_mat[pos_value])
    
    if (l == 1 && c == 1){
      data <- data.frame("id" = paste(as.character(l), as.character(c), sep = '.'), 'x' = x1, 'y' = y1)
    }
    for (i in 1:n){
      new_row <- data.frame("id" = paste(as.character(l), as.character(c), sep = '.'), 'x' = x1, 'y' = y1)
      data <- rbind(data, new_row)
    }
    end_loop_time <- Sys.time()
  }
  time_passed = end_loop_time - start_time
  expected_time = error * (end_loop_time - start_time) * ((total_lines - nrow(data))/nrow(data))
  message('Time passed: ', format(time_passed, format =  '%M:%S'), ' ... Expectation: ', 
          format(expected_time, format = '%M:%S'), '.')
}

# Save into a CSV. You don't want to lose it after all this work!!!
write.csv(data, 'map_data_new.csv')










