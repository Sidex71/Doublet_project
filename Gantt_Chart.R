###########################
library(scales)
library(ggplot2)
# Example data
tasks <- c("Task 1 \n yyy", "Task 2", "Task 3")
start_dates <- as.Date(c("2022-01-01", "2022-02-01", "2022-03-01"))
durations <- c(10, 15, 20)  # Duration of each task in days

# Create data frame
df <- data.frame(Task = tasks, StartDate = start_dates, Duration = durations)
df$EndDate <- df$StartDate + df$Duration
# Example data frame with additional column for colors
df$Color <- c("blue", "green", "red")  # Example colors for each task

# Create the Gantt chart with custom colors
ggplot(df, aes(x = StartDate, xend = EndDate, y = Task, yend = Task, color = Color)) +
  geom_segment(size = 10) +
  scale_x_date(date_breaks = "1 week", date_labels = "%b %d") +
  labs(x = "Date", y = "Task") +
  theme_minimal()

##################################
tasks <- c("Development of doublet detection tool \n (Thesis Chapter 2)", 
           "Cancer-Immune interaction modelling \n (Thesis Chapter 3)", 
           'Immune-Stromal interaction modelling \n (Thesis Chapter 3)', 
           'Research Placement',
           'Spatial transcriptomics validation \n (Thesis Chapter 4)', 
           'Thesis writing and submission')
start_dates <- as.Date(c("2024-03-01", "2024-06-01", 
                         "2025-01-01", '2025-07-01', '2025-09-01', 
                         '2024-03-01'))
durations <- c(90, 180, 180, 90, 180, 850)  # Duration of each task in days
df <- data.frame(Task = tasks, StartDate = start_dates, Duration = durations)
df$EndDate <- df$StartDate + df$Duration
# Example data frame with additional column for colors
df$Color <- c("blue", "green", "red", 'blue', 'green', 'red')  # Example colors for each task

# Reorder Task based on StartDate in descending order
df$Task <- factor(df$Task, levels = df$Task[order(df$EndDate, decreasing = TRUE)])

# Plot the Gantt chart with year labels on the x-axis
p1<-ggplot(df, aes(x = StartDate, xend = EndDate, y = Task, yend = Task, color = Color)) +
  geom_segment(size = 10) +
  scale_x_date(date_breaks = "2 month", date_labels = "%Y-%b") +
  labs(x = " ", y = "Tasks") +
  ggtitle('Gantt Chart showing doctoral research plans') +
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 45),
        axis.text.x = element_text(angle = 45, hjust = 1, size =40, face= "bold"),
        axis.text.y = element_text(size =40, face= "bold"), 
        plot.title = element_text(hjust = 0.5, size = 45, face = 'bold'),
        panel.grid.major = element_line(color = "black", size = 0.3),  # Set major grid lines
        panel.grid.minor = element_line(color = "black", size = 0.3))+  # Set minor grid lines) +  # Center title
  NoLegend()
p1
png("/mnt/Data/shameed/Gannt_Chart.png", width = 32.5, height = 16.5, units = 'in', res = 600)
p1
dev.off()

