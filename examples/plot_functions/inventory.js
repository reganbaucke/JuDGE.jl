<script src="https://cdn.plot.ly/plotly-2.1.0.min.js"></script>
<script>
function plotly_graph(data) {
  if (data==null) {
    document.getElementById("plotly").innerHTML="";
  }
  else {
    Plotly.newPlot('plotly', {
      data: data,
      layout: {
                width: 600,
                height: 400,
                margin: {
                  t:25
                },
                yaxis2: {
                  overlaying: 'y',
                  side: 'right'
                }
              }
    });
  }
}
</script>
