function(el, x, data) {
  // Each menu combination uses precomputed results from the public R function.
  // This is a finite teaching example, not a browser reimplementation of ocean3d.
  var P = window.Plotly;
  el.setAttribute('role', 'img');
  el.setAttribute('aria-label', 'Interactive example of ' + data.view + ' point selection; numerical results are also printed above the figure.');
  var controls = document.createElement('div');
  controls.style.cssText = 'display:flex;flex-wrap:wrap;gap:12px;margin:16px 0;';
  var status = document.createElement('div');
  status.setAttribute('aria-live', 'polite'); status.style.margin = '8px 0';
  var diagnostics = document.createElement('details');
  var summary = document.createElement('summary'); summary.textContent = 'Messages from this R extraction';
  var details = document.createElement('pre');
  details.style.cssText = 'white-space:pre-wrap;overflow-wrap:anywhere;';
  diagnostics.append(summary, details);
  el.parentNode.insertBefore(controls, el); el.parentNode.insertBefore(status, el);
  var caption = el.nextElementSibling;
  while (caption && caption.tagName === 'SCRIPT') caption = caption.nextElementSibling;
  if (caption && caption.tagName === 'P' && /^Figure\s+\d+\./.test(caption.textContent.trim())) {
    caption.classList.add('caption');
    caption.id = el.id + '-caption';
    el.setAttribute('aria-describedby', caption.id);
    caption.insertAdjacentElement('afterend', diagnostics);
  } else {
    el.parentNode.insertBefore(diagnostics, el.nextSibling);
  }
  var textColour = getComputedStyle(el).color, font = getComputedStyle(el).fontFamily;
  function f(n, digits) {return n === null || n === undefined ? 'NA' : Number(n).toFixed(digits === undefined ? 2 : digits);}
  function select(label, values) {
    var wrapper = document.createElement('label'); wrapper.className = 'form-label';
    wrapper.style.cssText = 'flex:1 1 145px;'; wrapper.append(document.createTextNode(label));
    var input = document.createElement('select'); input.className = 'form-select';
    values.forEach(function(v) {var option = document.createElement('option');
      option.value = v[0]; option.textContent = v[1]; input.append(option);});
    wrapper.append(input); controls.append(wrapper); return input;
  }
  function conditionText(record) {
    var lines = [].concat(record.warnings || [], record.messages || []);
    details.textContent = lines.length ? lines.join('\n\n') : 'No report or warning was needed for this selection.';
  }
  function layout() {return {height: 500, margin: {l: 65, r: 24, t: 20, b: 55},
    paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
    font: {family: font, color: textColour, size: 12},
    xaxis: {gridcolor: '#dddddd', zerolinecolor: '#cccccc'},
    yaxis: {gridcolor: '#dddddd', zerolinecolor: '#cccccc'}, showlegend: false};}
  function line(x, y, colour, width) {return {type: 'scatter', mode: 'lines', x: x, y: y,
    line: {color: colour, width: width || 1}, hoverinfo: 'skip', showlegend: false};}
  var draw, inputs = [], camera = {eye: {x: 1.8, y: -1.8, z: 1.1}};
  if (data.view === 'horizontal') {
    var location = select('Observation', data.locations.map(function(r) {return [r.key, r.label];}));
    var policy = select('Spatial matching', [['within', 'Within a cell'], ['nearest', 'Nearest cell centre']]);
    var limit = select('Nearest limit (km)', [['1', '1 km'], ['20', '20 km']]); limit.value = '20';
    inputs = [location, policy, limit];
    draw = function() {
      limit.disabled = policy.value === 'within';
      var record = data.horizontal.find(function(r) {return r.location === location.value && r.policy === policy.value && String(r.limit) === limit.value;});
      var selected = data.cells[record.cell - 1], traces = [];
      var values = data.cells.map(function(c) {return c.value;}), lo = Math.min.apply(null, values), hi = Math.max.apply(null, values);
      data.cells.forEach(function(c) {
        var p = (c.value - lo) / (hi - lo), stops = [[68,1,84],[59,82,139],[33,145,140],[94,201,98],[253,231,37]];
        var pos = p * 4, a = Math.floor(pos), b = Math.min(a + 1, 4);
        var rgb = stops[a].map(function(v, i) {return Math.round(v + (stops[b][i] - v) * (pos - a));});
        traces.push({type: 'scatter', mode: 'lines', x: c.polygon.map(function(p) {return p[0];}),
          y: c.polygon.map(function(p) {return p[1];}), fill: 'toself', fillcolor: 'rgb(' + rgb.join(',') + ')',
          line: {color: '#bbbbbb', width: 0.6}, hoveron: 'fills',
          text: 'Cell ' + c.id + ': ' + f(c.value) + ' degrees C<br>Centre: ' + f(c.lon, 5) + ', ' + f(c.lat, 5),
          hovertemplate: '%{text}<extra></extra>'});
      });
      if (record.value !== null) {
        traces.push(line(selected.polygon.map(function(p) {return p[0];}), selected.polygon.map(function(p) {return p[1];}), textColour, 3));
        traces.push(line([record.point[0], selected.centre[0]], [record.point[1], selected.centre[1]], textColour, 2));
        traces.push({type: 'scatter', mode: 'markers', x: [selected.centre[0]], y: [selected.centre[1]],
          marker: {size: 8, color: textColour}, hovertemplate: 'Selected cell centre<extra></extra>'});
      }
      traces.push({type: 'scatter', mode: 'markers', x: [record.point[0]], y: [record.point[1]],
        marker: {size: 14, symbol: 'x', color: 'red'}, hovertemplate: 'Requested observation<extra></extra>'});
      status.textContent = 'Returned temperature: ' + f(record.value) + ' degrees C | distance to nearest centre: ' + f(record.distance_km) + ' km';
      conditionText(record);
      var l = layout(); l.xaxis.title = {text: 'East of original observation (km)'};
      l.yaxis.title = {text: 'North of original observation (km)'};
      l.xaxis.range = [-23, 27]; l.yaxis.range = [-23, 23];
      l.xaxis.constrain = 'domain'; l.yaxis.constrain = 'domain'; l.yaxis.scaleanchor = 'x';
      return P.react(el, traces, l, {responsive: true, displaylogo: false});
    };
  } else if (data.view === 'vertical') {
    var date = select('Date (UTC)', data.profiles.map(function(p) {return [p.date, p.date];})); date.value = '2026-09-20';
    var target = select('Requested depth', [['60', '60 m'], ['100', '100 m'], ['800', '800 m (outside depth range)']]);
    var method = select('Depth method', [['nearest', 'Nearest valid layer'], ['surface', 'Shallowest layer'], ['bottom', 'Deepest valid layer'], ['all', 'All three selections']]);
    var view = select('View', [['profile', 'Temperature profile'], ['column', '3D selected cell column']]);
    inputs = [date, target, method, view];
    draw = function() {
      var record = data.vertical.find(function(r) {return r.date === date.value && String(r.target) === target.value;});
      var profile = data.profiles.find(function(p) {return p.date === date.value;});
      var methods = method.value === 'all' ? ['nearest', 'surface', 'bottom'] : [method.value];
      var indices = methods.map(function(m) {return record[m] - 1;});
      conditionText(record);
      // The stored reports come from method = all, so all three outputs are named.
      summary.textContent = 'Messages from the R extraction with method = "all"';
      status.textContent = methods.map(function(m, i) {return m + ': ' + f(data.depths[indices[i]]) + ' m, ' + f(profile.values[indices[i]]) + ' degrees C';}).join(' | ');
      var traces = [], l = layout();
      if (view.value === 'profile') {
        traces.push({type: 'scatter', mode: 'lines+markers', x: profile.values, y: data.depths,
          marker: {size: 5, color: '#777777'}, line: {color: '#777777', width: 1},
          hovertemplate: '%{y:.2f} m: %{x:.2f} degrees C<extra>Original layer</extra>'});
        var range = [Math.min.apply(null, profile.values) - 0.5, Math.max.apply(null, profile.values) + 0.5];
        traces.push(line(range, [record.target, record.target], 'red', 2));
        var colours = ['darkorange', 'steelblue', 'purple'];
        indices.forEach(function(index, i) {traces.push({type: 'scatter', mode: 'markers+text',
          x: [profile.values[index]], y: [data.depths[index]], text: [methods[i]], textposition: methods[i] === 'surface' ? 'bottom left' : 'top right',
          marker: {size: 11, color: colours[i], symbol: ['circle', 'diamond', 'square'][i]},
          hovertemplate: methods[i] + ': %{y:.2f} m, %{x:.2f} degrees C<extra></extra>'});});
        l.xaxis.title = {text: 'Potential temperature (degrees C)'};
        l.yaxis.title = {text: 'Depth (m, positive down)'};
        l.yaxis.range = [Math.max(430, record.target + 30), -15];
      } else {
        var cell = data.cell, polygon = cell.polygon.slice(0, -1), n = polygon.length;
        var xx = [], yy = [], zz = [], ii = [], jj = [], kk = [], intensity = [];
        data.depths.forEach(function(z, k) {
          var start = xx.length;
          data.depth_bounds[k].forEach(function(bound) {polygon.forEach(function(p) {
            xx.push(p[0]); yy.push(p[1]); zz.push(bound); intensity.push(profile.values[k]);
          });});
          for (var q = 1; q < n - 1; q++) {ii.push(start, start+n); jj.push(start+q, start+n+q+1); kk.push(start+q+1, start+n+q);}
          for (var q = 0; q < n; q++) {var next=(q+1)%n;ii.push(start+q,start+q);jj.push(start+next,start+n+next);kk.push(start+n+next,start+n+q);}
        });
        var all = data.profiles.flatMap(function(p) {return p.values;});
        traces.push({type: 'mesh3d', x: xx, y: yy, z: zz, i: ii, j: jj, k: kk,
          intensity: intensity, colorscale: 'Viridis', cmin: Math.min.apply(null, all), cmax: Math.max.apply(null, all),
          flatshading: true, opacity: 0.85, colorbar: {title: {text: 'degrees C'}, thickness: 12, len: 0.7},
          hovertemplate: '%{z:.2f} m (display support)<extra></extra>', lighting: {ambient: 0.8, specular: 0}});
        indices.forEach(function(index, i) {traces.push({type: 'scatter3d', mode: 'markers',
          x: [cell.centre[0]], y: [cell.centre[1]], z: [data.depths[index]], text: [methods[i]],
          marker: {size: 7, color: 'darkorange'}, hovertemplate: methods[i] + ': %{z:.2f} m<extra></extra>'});});
        traces.push({type: 'scatter3d', mode: 'markers', x: [0], y: [0], z: [record.target],
          marker: {size: 6, color: 'red'}, hovertemplate: 'Requested depth: %{z} m<extra></extra>'});
        l.margin = {l: 0, r: 10, t: 5, b: 0};
        var selectedAnnotations = indices.map(function(index, i) {return {
          x:cell.centre[0],y:cell.centre[1],z:data.depths[index],
          text:methods[i]+': '+f(data.depths[index])+' m',showarrow:true,
          ax:45,ay:methods[i]==='bottom'?20:-15,bgcolor:'rgba(255,255,255,0.9)',
          font:{color:textColour},arrowcolor:textColour
        };});
        l.scene = {camera: camera, aspectmode: 'manual', aspectratio: {x: 0.65, y: 0.65, z: 1.5},
          xaxis: {title: {text: 'East (km)'}}, yaxis: {title: {text: 'North (km)'}},
          zaxis: {title: {text: 'Depth (m)'}, range: [Math.max(430, record.target+30), 0]},
          annotations: selectedAnnotations.concat([{x:0,y:0,z:record.target,text:'Requested: '+record.target+' m',
            showarrow:true,ax:-45,ay:-45,bgcolor:'rgba(255,255,255,0.9)',font:{color:textColour},arrowcolor:'red'}])};
        l.uirevision = 'point-column';
      }
      return P.react(el, traces, l, {responsive: true, displaylogo: false});
    };
  } else {
    var requests = Array.from(new Set(data.temporal.map(function(r) {return r.request;})));
    var request = select('Requested time (UTC)', requests.map(function(r) {return [r,r];})); request.value = requests[1];
    var precision = select('Comparison precision', [['datetime','Datetime'],['daily','Daily labels']]);
    var matching = select('Time matching', [['nearest','Nearest'],['exact','Exact']]);
    var tolerance = select('Maximum gap (days)', [['none','No limit'],['0.1','0.1 days (2.4 h)'],['1','1 day']]);
    inputs = [request,precision,matching,tolerance];
    draw = function() {
      var record = data.temporal.find(function(r) {return r.request===request.value&&r.precision===precision.value&&r.matching===matching.value&&r.limit===tolerance.value;});
      conditionText(record);
      status.textContent = 'Returned temperature: '+f(record.value)+' degrees C | selected record: '+(record.selected || 'none')+(record.selected?' 00:00 UTC':'');
      var available = data.profiles.map(function(p) {return p.date+'T00:00:00Z';});
      var comparison = record.precision === 'daily' ? record.request.slice(0,10)+'T00:00:00Z' : record.request.replace(' ','T')+'Z';
      var traces=[{type:'scatter',mode:'markers+text',x:available,y:available.map(function(){return 1;}),
        text:data.profiles.map(function(p){return p.date.slice(5);}),textposition:'top center',marker:{size:12,color:'steelblue'},
        hovertemplate:'Available: %{x}<extra></extra>'},
        {type:'scatter',mode:'markers',x:[record.request.replace(' ','T')+'Z'],y:[0.45],
          marker:{size:15,color:'red',symbol:'x'},hovertemplate:'Requested instant: %{x}<extra></extra>'},
        {type:'scatter',mode:'markers',x:[comparison],y:[0.2],
          marker:{size:12,color:'darkorange',symbol:'diamond'},hovertemplate:'Comparison value: %{x}<extra></extra>'}];
      if(record.selected) traces.push({type:'scatter',mode:'markers',x:[record.selected+'T00:00:00Z'],y:[1],
        marker:{size:22,color:'rgba(0,0,0,0)',line:{color:textColour,width:3}},hovertemplate:'Selected record: %{x}<extra></extra>'});
      var l=layout();l.height=360;l.margin={l:95,r:25,t:30,b:65};
      l.xaxis={type:'date',title:{text:'UTC date and time'},tickformat:'%b %d<br>%H:%M',nticks:5,gridcolor:'#dddddd',
        range:['2026-09-18T12:00:00Z',record.request.slice(0,10)==='2026-09-25'?'2026-09-25T12:00:00Z':'2026-09-21T12:00:00Z']};
      l.yaxis={range:[0,1.3],tickvals:[0.2,0.45,1],ticktext:['Compared','Requested','Available'],showgrid:false,zeroline:false};
      return P.react(el,traces,l,{responsive:true,displaylogo:false});
    };
  }
  var busy=false,queued=false;
  async function render(){if(busy){queued=true;return;}busy=true;try{await draw();}finally{busy=false;if(queued){queued=false;render();}}}
  inputs.forEach(function(input){input.addEventListener('change',render);});
  el.on('plotly_relayout',function(event){if(event['scene.camera'])camera=event['scene.camera'];});
  render();
}
