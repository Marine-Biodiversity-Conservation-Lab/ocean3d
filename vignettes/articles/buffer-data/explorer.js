function(el, x, data) {
  // The data and all Plotly dependencies are included in the built article.
  // Controls change the display only; they do not rerun an extraction.
  var P = window.Plotly;
  var controls = document.createElement('div');
  controls.style.cssText = 'display:flex;flex-wrap:wrap;gap:12px;margin:16px 0;';
  var status = document.createElement('div');
  status.setAttribute('aria-live', 'polite');
  status.style.margin = '8px 0';
  el.parentNode.insertBefore(controls, el);
  el.parentNode.insertBefore(status, el);
  function select(label, values) {
    var wrapper = document.createElement('label');
    wrapper.className = 'form-label';
    wrapper.style.cssText = 'flex:1 1 150px;';
    wrapper.append(document.createTextNode(label));
    var input = document.createElement('select');
    input.className = 'form-select';
    values.forEach(function(v) {
      var option = document.createElement('option');
      option.value = v[0]; option.textContent = v[1]; input.append(option);
    });
    wrapper.append(input); controls.append(wrapper); return input;
  }
  function f(value, digits) { return Number(value).toFixed(digits === undefined ? 2 : digits); }
  var view = select('View', [['3d', '3D buffer'], ['plan', 'Horizontal footprint']]);
  var time = select('Date (UTC)', data.times.map(function(t) {return [t.index, t.label];}));
  time.value = data.times[1].index;
  var layer = select('Depth layer', [['all', 'All layers: 50-70 m']].concat(
    data.depths.map(function(d) {return [d.index, f(d.centre) + ' m'];})));
  var colour = select('Colour', [['temperature', 'Potential temperature'], ['weight', 'Contribution weight']]);
  var weighting = select('Weighting', [['weighted', 'Weighted'], ['equal', 'Equal']]);
  var allValues = data.samples.map(function(s) {return s.value;});
  var minimum = Math.min.apply(null, allValues), maximum = Math.max.apply(null, allValues);
  var maxWeight = Math.max.apply(null, data.samples.map(function(s) {return s.weight * 100;}));
  var camera = {eye: {x: 1.6, y: -1.8, z: 1.35}};
  var font = getComputedStyle(el).fontFamily;
  var textColour = getComputedStyle(el).color;
  var busy = false, queued = false;
  function depth(index) { return data.depths.find(function(d) {return d.index === index;}); }
  function weight(s) {return weighting.value === 'equal' ? 1 / data.samples.length : s.weight;}
  function value(s) {return colour.value === 'temperature' ? s.value : 100 * weight(s);}
  function description(s, cell) {
    var d = depth(s.depth_index);
    return 'Cell ' + cell.id + ' | ' + f(s.value) + ' degrees C<br>' +
      'Layer centre: ' + f(d.centre) + ' m; used: ' + f(d.lower) + '-' + f(d.upper) + ' m<br>' +
      'Weight: ' + f(100 * weight(s), 4) + '% of all 450 values';
  }
  function axis(title, range) {
    return {title: {text: title}, range: range, gridcolor: '#dddddd',
      zerolinecolor: '#cccccc', showbackground: false, tickfont: {color: textColour}};
  }
  async function render() {
    if (busy) {queued = true; return;}
    busy = true;
    try {
      var samples = data.samples.filter(function(s) {
        return String(s.time_index) === time.value &&
          (layer.value === 'all' || String(s.depth_index) === layer.value);
      });
      var traces = [], isTemp = colour.value === 'temperature';
      var scale = isTemp ? 'Viridis' : 'Blues';
      var cmin = isTemp ? minimum : 0, cmax = isTemp ? maximum : maxWeight;
      var date = data.times.find(function(t) {return String(t.index) === time.value;});
      status.textContent = 'Mean across all 450 values: ' +
        f(weighting.value === 'weighted' ? data.weighted : data.equal, 5) +
        ' degrees C | showing ' + samples.length + ' values on ' + date.label.slice(0, 10) + ' UTC';
      var layout = {height: 560, margin: {l: 0, r: 20, t: 10, b: 0},
        paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
        font: {family: font, color: textColour, size: 12}, showlegend: false};
      if (view.value === '3d') {
        var xx = [], yy = [], zz = [], ii = [], jj = [], kk = [], intensity = [], text = [];
        samples.forEach(function(s) {
          var cell = data.cells[s.cell - 1], d = depth(s.depth_index);
          var ring = cell.polygon.slice(0, -1), n = ring.length, start = xx.length;
          [d.lower, d.upper].forEach(function(z) {ring.forEach(function(p) {
            xx.push(p[0]); yy.push(p[1]); zz.push(z);
            intensity.push(value(s)); text.push(description(s, cell));
          });});
          for (var q = 1; q < n - 1; q++) {
            ii.push(start, start + n); jj.push(start + q, start + n + q + 1);
            kk.push(start + q + 1, start + n + q);
          }
          for (var q = 0; q < n; q++) {
            var next = (q + 1) % n;
            ii.push(start + q, start + q); jj.push(start + next, start + n + next);
            kk.push(start + n + next, start + n + q);
          }
        });
        traces.push({type: 'mesh3d', x: xx, y: yy, z: zz, i: ii, j: jj, k: kk,
          intensity: intensity, colorscale: scale, cmin: cmin, cmax: cmax,
          colorbar: {title: {text: isTemp ? 'degrees C' : 'Weight (%)'}, thickness: 12, len: 0.65},
          flatshading: true, opacity: 0.94, text: text,
          hovertemplate: '%{text}<extra></extra>', lighting: {ambient: 0.8, diffuse: 0.5, specular: 0}});
        [50, 70].forEach(function(z) {
          var angles = Array.from({length: 121}, function(_, i) {return 2 * Math.PI * i / 120;});
          traces.push({type: 'scatter3d', mode: 'lines',
            x: angles.map(function(a) {return 10 * Math.cos(a);}),
            y: angles.map(function(a) {return 10 * Math.sin(a);}), z: angles.map(function() {return z;}),
            line: {color: textColour, width: 3}, hoverinfo: 'skip'});
        });
        traces.push({type: 'scatter3d', mode: 'markers', x: [0], y: [0], z: [60],
          marker: {color: 'red', size: 6, symbol: 'diamond'},
          hovertemplate: 'Observation: 0 E, 38 N, 60 m<extra></extra>'});
        layout.scene = {xaxis: axis('East (km)', [-12, 12]), yaxis: axis('North (km)', [-12, 12]),
          zaxis: axis('Depth (m)', [72, 48]), camera: camera,
          aspectmode: 'manual', aspectratio: {x: 1, y: 1, z: 0.65},
          annotations: [{x: 0, y: 0, z: 60, text: 'Observation: 60 m', showarrow: true,
            ax: 45, ay: -55, arrowcolor: textColour, font: {color: textColour}}]};
        layout.uirevision = 'buffer-3d';
      } else {
        data.cells.forEach(function(cell) {
          var selected = samples.filter(function(s) {return s.cell === cell.id;});
          if (!selected.length) return;
          var total = selected.reduce(function(a, s) {return a + weight(s);}, 0);
          var temp = selected.reduce(function(a, s) {return a + s.value * weight(s);}, 0) / total;
          // Filled triangles give the clipped footprint a constant value per cell.
          var ring = cell.polygon;
          var display = isTemp ? temp : 100 * total;
          var lo = isTemp ? minimum : 0;
          var hi = isTemp ? maximum : maxWeight * data.depths.length;
          var proportion = Math.max(0, Math.min(1, (display - lo) / (hi - lo)));
          var palette = isTemp ? [[68,1,84],[59,82,139],[33,145,140],[94,201,98],[253,231,37]] : [[247,251,255],[107,174,214],[8,48,107]];
          var pos = proportion * (palette.length - 1), a = Math.floor(pos), b = Math.min(a + 1, palette.length - 1);
          var rgb = palette[a].map(function(v, i) {return Math.round(v + (palette[b][i] - v) * (pos - a));});
          traces.push({type: 'scatter', mode: 'lines', x: cell.full.map(function(p) {return p[0];}),
            y: cell.full.map(function(p) {return p[1];}), line: {color: '#bbbbbb', width: 1}, hoverinfo: 'skip'});
          traces.push({type: 'scatter', mode: 'lines', x: ring.map(function(p) {return p[0];}),
            y: ring.map(function(p) {return p[1];}), fill: 'toself', fillcolor: 'rgb(' + rgb.join(',') + ')',
            line: {color: '#777777', width: 0.5}, hoveron: 'fills',
            text: 'Cell ' + cell.id + '<br>' + (isTemp ? f(display) + ' degrees C' : f(display, 4) + '% of total weight'),
            hovertemplate: '%{text}<extra></extra>'});
        });
        traces.push({type: 'scatter', mode: 'markers', x: [0], y: [0],
          marker: {color: 'red', size: 12, symbol: 'x'}, hovertemplate: 'Observation<extra></extra>'});
        layout.xaxis = axis('East of observation (km)', [-13, 13]);
        layout.yaxis = Object.assign(axis('North of observation (km)', [-13, 13]), {scaleanchor: 'x', scaleratio: 1});
        layout.margin = {l: 60, r: 20, t: 10, b: 50}; layout.uirevision = 'buffer-plan';
      }
      await P.react(el, traces, layout, {responsive: true, displaylogo: false});
    } finally {busy = false; if (queued) {queued = false; render();}}
  }
  [view, time, layer, colour, weighting].forEach(function(input) {input.addEventListener('change', render);});
  el.on('plotly_relayout', function(event) {if (event['scene.camera']) camera = event['scene.camera'];});
  render();
}
