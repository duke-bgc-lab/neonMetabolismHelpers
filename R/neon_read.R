read_all_neon_feathers <- function(file_path, by_site = TRUE){

    if(by_site == TRUE){


        # full paht ot each NEON site directory in data location
        site_dirs <- list.files(file_path, full.names = TRUE)
        # get full path of all files inside all NEON site directories
        neon_files <- list.files(site_dirs, full.names = TRUE)
        # get just file names fpr each of these files
        file_names <- list.files(site_dirs)

        neon_list <- lapply(neon_files, read_feather)
        names(neon_list) <- file_names

        return(neon_list)

    } else{

        neon_files <- list.files(file_path, full.names = TRUE, recursive = TRUE)
        file_names <- str_split_fixed(neon_files, '.feather', n = Inf)[,1]
        file_names <- str_split_fixed(file_names, '[/]', n = Inf)
        inv_file_names <- file_names[,dim(file_names)[2]]
        site_name <- file_names[,dim(file_names)[2]-1]

        inv_file_names <- unique(inv_file_names)
        site_names <- unique(site_name)

        #categoricalCodes
        #readme
        #validation
        #variables

        data_files <-  inv_file_names[!grepl('categoricalCodes|readme|validation|variables', inv_file_names)]

        final_list <- list()
        for(i in 1:length(data_files)){

            all_files <- tibble()
            for(s in 1:length(site_names)){

                path <- glue('{f}/{s}/{i}.feather',
                             f = file_path,
                             i = data_files[i],
                             s = site_names[s])

                one_site <- read_feather(path)

                all_files <- rbind(all_files, one_site)
            }

            all_site_files <- list(all_files)

            names(all_site_files) <- data_files[i]

            final_list <- c(final_list, all_site_files)
        }

        meta_data_files <-  inv_file_names[grepl('categoricalCodes|readme|validation|variables', inv_file_names)]

        for(i in 1:length(meta_data_files)){

            path <- glue('{f}/{s}/{i}.feather',
                         f = file_path,
                         i = meta_data_files[i],
                         s = site_names[1])

            meta_file <- read_feather(path)

            meta_list <- list(meta_file)
            names(meta_list) <- meta_data_files[i]

            final_list <- c(final_list, meta_list)
        }
        return(final_list)
    }

}

read_neon_feathers <- function(file_path, by_site){

    if(by_site == TRUE){

        neon_files <- list.files(file_path, full.names = TRUE)

        file_names <- list.files(file_path)

        neon_list <- lapply(neon_files, read_feather)

        names(neon_list) <- file_names

        return(neon_list)

    } else{

        neon_files <- list.files(file_path, full.names = TRUE, recursive = TRUE)
        file_names <- str_split_fixed(neon_files, '.feather', n = Inf)[,1]
        file_names <- str_split_fixed(file_names, '[/]', n = Inf)
        inv_file_names <- file_names[,dim(file_names)[2]]
        site_name <- file_names[,dim(file_names)[2]-1]

        inv_file_names <- unique(inv_file_names)
        site_names <- unique(site_name)

        #categoricalCodes
        #readme
        #validation
        #variables

        data_files <-  inv_file_names[!grepl('categoricalCodes|readme|validation|variables', inv_file_names)]

        final_list <- list()
        for(i in 1:length(data_files)){

            all_files <- tibble()
            for(s in 1:length(site_names)){

                path <- glue('{f}/{s}/{i}.feather',
                             f = file_path,
                             i = data_files[i],
                             s = site_names[s])

                one_site <- read_feather(path)

                all_files <- rbind(all_files, one_site)
            }

            all_site_files <- list(all_files)

            names(all_site_files) <- data_files[i]

            final_list <- c(final_list, all_site_files)
        }

        meta_data_files <-  inv_file_names[grepl('categoricalCodes|readme|validation|variables', inv_file_names)]

        for(i in 1:length(meta_data_files)){

            path <- glue('{f}/{s}/{i}.feather',
                         f = file_path,
                         i = meta_data_files[i],
                         s = site_names[1])

            meta_file <- read_feather(path)

            meta_list <- list(meta_file)
            names(meta_list) <- meta_data_files[i]

            final_list <- c(final_list, meta_list)
        }
        return(final_list)
    }

}
