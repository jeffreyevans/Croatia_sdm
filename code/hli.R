hli <- function(aspect, slope, latitude, direct = FALSE, 
                units = c("degrees", "radians"), 
				hemisphere = c("northern", "southern"),
				force.hemisphere = TRUE, equation = 1){
  if (units[1] == "degrees") {
    message("Converting degrees to radians")
      aspect <- aspect / 180 * pi
        slope <- slope / 180 * pi
      aspect[slope == 0] <- 0
	if(missing(latitude)) {
	  e <- st_as_sf(as.polygons(ext(ref)))
	    st_crs(e) <- st_crs(ref)
	  latitude <- as.numeric(sf::st_coordinates(sf::st_centroid(
	                         st_transform(e, st_crs(4326))))[,2])
	}
    latitude <- latitude / 180 * pi
  } else if(units[1] == "radians") {
    message("Assuming data is in radians")
  }  
      
   if(hemisphere[1] == "northern" & force.hemisphere == TRUE){ 
	  message("Using folded aspect equation for Northern hemisphere") 	
	    # Folded Aspect Northern Hemisphere  (180 - (Aspect – 225) )
		if(!direct) {
          A <- abs(pi - abs (aspect - (5 * pi / 4))) 
        } else {
		  A <- pi - abs (aspect - pi)
		}
	} else if(hemisphere[1] == "southern" & force.hemisphere == TRUE) {
	  message("Using folded aspect equation for Southern hemisphere") 		
	  # Folded Aspect Southern Hemisphere  ( 180 - ( Aspect – 315) )  
	  A <- abs(3.141593 - abs(aspect - 5.497791))
	}	  
  S <- slope
  L <- if (length (latitude) == 1) rep (latitude, length (A)) else latitude 
  if (equation == 1) {
    res <- exp(-1.467 + 1.582 * cos(L) * cos(S) -1.500 * cos(A) * sin(S) * 
	           sin(L) -0.262 * sin(L) * sin(S) + 0.607 * sin(A) * sin(S))
  } else if(equation == 2) { 
    res <- exp(-1.236 + 1.350 * cos(L) * cos(S) - 1.376 * cos(A) * sin(S) * sin(L) - 
               0.331 * sin(L) * sin(S) + 0.375 * sin(A) * sin(S))
  } else if (equation == 3) { 
    res <- 0.339 + (0.808 * cos(L) * cos(S))  - (0.196 * sin(L) * sin(S))  - 
		   (0.482 * cos(A) * sin(S))
  }
  return (res)
}
