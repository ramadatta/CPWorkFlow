https://www.biostars.org/p/9490989/

### Adding multiple timepoint on dumbbell plot


```
read.table(text="whowid  gasample1_weeks gasample2_weeks gaoutcome_weeks
36790   19  24  40
36821   19  32  39
36912   19  32  38
36916   19  26  37
36941   19  24  39
36994   18  28  44
37043   19  24  37
37046   16  24  38
37063   19  32  40
37068   19  25  40
37117   19  24  42
37120   16  24  39
37122   16  32  37
37168   18  28  38", header = TRUE) -> xdf

xdf$whowid <- factor(xdf$whowid)

ggplot(data = xdf) +
  geom_dumbbell(
    aes(y = whowid, x=gasample1_weeks, xend=gaoutcome_weeks),
    size_x = 3, size_xend = 3
  ) +
  geom_point(
    aes(y = whowid, x = gasample2_weeks),
    size = 3, color = "steelblue"
  )
  ```
