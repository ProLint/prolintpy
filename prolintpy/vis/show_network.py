def show_network(filename, width=800, height=600):
    """Visualize lipid-protein contacts using a network graph.

    Parameters
    ----------

    filename : str
        The name of the json file in which network information is stored. This
        is the output of pl.network() call.

    width : int
        The width of the window

    height : int
        The height of the window.

    """

    from IPython.core.display import HTML
    from string import Template, ascii_uppercase, ascii_lowercase, digits
    import json
    import random
    import string

    rand_str = ''.join(random.choices(ascii_uppercase + digits + ascii_lowercase, k=10))

    with open(filename) as f:
        network = json.load(f)

    HTML('<script src="./d3.min.js"></script>')

    js_text_template = Template('''

    var clientHeight = document.getElementById('$network').clientHeight;
    var clientWidth = document.getElementById('$network').clientWidth;

    var root = $python_data ;

    var width = clientWidth,
        height = clientHeight,
        root,
        radius = 6;

    function loadDATA(gpcr) {
        var force = d3.layout.force()
            .size([width, height])
            .on("tick", tick);


        force.charge(function (d) {
            if (d.id <= root.graph.lip_end) {
                return -500
            }
            return -5000
        });


        force.linkDistance(function (d) {
            if ((d.source.id <= root.graph.lip_end) && (d.target.id >= root.graph.res_start + 1)) {
                return 300
            }
            return 180
        });


        var svg = d3.select("#$network").append("svg")
            .attr("width", width)
            .attr("height", height);


        var link = svg.selectAll(".link"),
            node = svg.selectAll(".node");

    function call_root(gpcr) {

        // root = gpcr;

        var a = width / 2
        var b = height / 2


        if (root.graph.lip_end <= 8) {

            var R = 50
            posX = [a, a + R, a - R, a, a, a + (R / 1.41), a + (R / 1.41), a - (R / 1.41), a - (R / 1.41)]
            posY = [b, b, b, b + R, b - R, b + (R / 1.41), b - (R / 1.41), b + (R / 1.41), b - (R / 1.41)]


        } else {

            var R = [50, 70, 90, 110, 30, 130, 140, 160, 180, 200];
            var new_r = 0
            var cp = 1
            var it = 1
            // Hendecagon alignment of lipid nodes
            angle = 0.57
            // angle = 0.785
            var posX = [a];
            var posY = [b];

            while (cp <= root.graph.lip_end) {

                var x = Math.cos(angle*it) * R[new_r];
                var y = Math.sin(angle*it) * R[new_r];
                posX.push(a + x)
                posY.push(b + y)

                if (it < 11) {
                    it++
                } else {
                    var it = 1
                    new_r++
                }
                cp++
            }
    }

        for (var i = 0; i < root.nodes.length; i++) {
            var node = root.nodes[i];

            if (node.id <= root.graph.lip_end) {
                node.x = posX[node.id]
                node.y = posY[node.id]
                node.fixed = true;
            }

            node.id = i;
            node.collapsing = 0;
            node.collapsed = false;
        }

        //Give links ids and initialize variables
        for (var i = 0; i < root.links.length; i++) {
            var link = root.links[i];
            link.source = root.nodes[link.source];
            link.target = root.nodes[link.target];
            link.id = i;
        }

        function hideChildren(d) {
            //check if link is from this node, and if so, collapse
            root.links.forEach(function (l) {
                if ((l.source.id >= root.graph.lip_end + 1) && (l.source.id <= root.graph.res_start)) {
                    l.target.collapsing = 1
                    l.source.collapsed = true
                }
            });

            update();

        }

        for (var i = 0; i < root.nodes.length; i++) {
            hideChildren(node)
        }

        update();

   // });
    };
    call_root(gpcr)

        function update() {

            d3.selectAll('text').remove()

            //Keep only the visible nodes
            var nodes = root.nodes.filter(function (d) {
                return d.collapsing == 0;
            });
            var links = root.links;

            //Keep only the visible links
            links = root.links.filter(function (d) {
                return d.source.collapsing == 0 && d.target.collapsing == 0;
            });


            force
                .nodes(nodes)
                .links(links)
                .friction(0.4)
                .start();

            // Update the links…
            link = link.data(links, function (d) {
                return d.id;
            });

            // Exit any old links.
            link.exit().remove();

            // Enter any new links.
            link.enter().insert("line", ".node")
                .attr("class", "link")
                .attr("x1", function (d) {
                    return d.source.x;
                })
                .attr("y1", function (d) {
                    return d.source.y;
                })
                .attr("x2", function (d) {
                    return d.target.x;
                })
                .attr("y2", function (d) {
                    return d.target.y;
                })
                // .style('stroke-width', 5)
                .style('stroke-width', function (d) {
                    return d.weight * 10;
                })
                .style("stroke-opacity", 0.6)
                .attr("stroke", "#b7b7b7")
                .attr("marker-end", "url(#arrow)");

            // Update the nodes…
            node = node.data(nodes, function (d) {
                return d.id;

            }).style("fill", color);

            // Exit any old nodes.
            node.exit().remove();

            // Enter any new nodes.
            node.enter().append("circle")
                .attr("class", "node")
                .attr("cx", function (d) {
                    return d.x;
                })
                .attr("cy", function (d) {
                    return d.y;
                })
                .attr("r", function (d) {
                    return d.values;
                })
                .style("fill", color)
                .style("cursor", "pointer")
                .on("click", click)
                .on("mouseover", mouseOver(.0))
                .on("mouseout", mouseOut)
                .call(force.drag);

            label = svg.selectAll(null)
                .data(nodes)
                .enter()
                .append("text")
                .text(function (d) {
                    return d.name;
                })
                .style("text-anchor", "middle")
                // .style("fill", "black")
                .style("fill", function(d) {
                    if (d.id <= root.graph.lip_end) {
                        return "#a03271"
                    }
                    return "black"
                })
                .style("font-family", "Arial")
                .style("font-weight", function(d) {
                    if (d.id <= root.graph.lip_end) {
                        return "bold"
                    }
                })
                .style("font-size", function (d) {
                    if (d.id <= root.graph.lip_end) {
                        return 10
                    }
                    return 12
                });

            // fade nodes on hover
            function mouseOver(opacity) {
                return function (d) {

                    var linkedByIndex = {};
                    root.links.forEach(function (d) {
                        linkedByIndex[d.source.id + "," + d.target.id] = 1;
                    });

                    function isConnected(a, b) {
                        return linkedByIndex[a.id + "," + b.id] || linkedByIndex[b.id + "," + a.id] || a.id == b.id;
                    }

                    // check all other nodes to see if they're connected
                    // to this one. if so, keep the opacity at 1, otherwise
                    // fade. Tweak text transparency/opacity
                    label.style("opacity", function (o) {
                        thisOpacity = isConnected(d, o) ? 1 : opacity;
                        return thisOpacity;
                    })

                    node.style("stroke-opacity", function (o) {
                        thisOpacity = isConnected(d, o) ? 0 : opacity;
                        return thisOpacity;
                    });
                    node.style("fill-opacity", function (o) {
                        thisOpacity = isConnected(d, o) ? 1 : opacity;
                        return thisOpacity;
                    });

                    // also style link accordingly
                    link.style("stroke-opacity", function (o) {
                        o.source.colour = "black"
                        return o.source === d || o.target === d ? 0.7 : opacity;
                    });

                    link.style("stroke", function (o) {
                        return o.source === d || o.target === d ? o.source.colour : "#b7b7b7";
                    });
                };
            }

            function mouseOut() {
                node.style("stroke-opacity", 1);
                node.style("fill-opacity", 1);
                link.style("stroke-opacity", 0.5);
                link.style("stroke", "#b7b7b7");
                label.style("opacity", 1);
            }
            // new new
            // if(abc) return;
        }

        function tick() {

            link.attr("x1", function (d) {
                    return d.source.x;
                })
                .attr("y1", function (d) {
                    return d.source.y;
                })
                .attr("x2", function (d) {
                    return d.target.x;
                })
                .attr("y2", function (d) {
                    return d.target.y;
                });

            // Bounding box example: requires radius
            node.attr("cx", function (d) {
                    return d.x = Math.max(radius, Math.min(width - radius, d.x));
                })
                .attr("cy", function (d) {
                    return d.y = Math.max(radius, Math.min(height - radius, d.y));
                });

            node
                .attr("cx", function (d) {
                    return d.x;
                })
                .attr("cy", function (d) {
                    return d.y;
                });

            label
                .attr("x", function (d) {
                    return d.x;
                })
                .attr("y", function (d) {
                    return d.y - 15;
                })

        }

        // Color leaf nodes light blue, and packages grey or black.
        function color(d) {
            if (d.id <= root.graph.lip_end) {
                // return d.collapsed = "#a03271"
                return d.collapsed ? "#d6c8ff" : d.children ? "#c6dbef" : "#a03271";
            }
            return d.collapsed ? "#4b4c7a" : d.children ? "#c6dbef" : "#d6c8ff";
        }

        // Toggle children on click.
        function click(d) {

                //check if link is from this node, and if so, collapse
                root.links.forEach(function (l) {

                    if (d.id <= root.graph.lip_end) {

                        if (d.collapsed) {
                            l.target.collapsing--;
                        } else {
                            l.target.collapsing++;
                        }
                    } else if (l.source.id == d.id) {

                        if (d.collapsed) {

                            l.target.collapsing--;

                        } else {

                            l.target.collapsing++;

                        }
                    }
                });

                d.collapsed = !d.collapsed;

            update();
        }
    }

    function loadGPCR() {
        loadDATA('')
    }
    loadGPCR()
    '''#.format(width, height)
    )

    html_template = Template(
    '<div id=' + rand_str + f' style="width:{width}px; height:{height}px"></div>'
    '''
    <script > $js_text </script>
    ''')

    js_text = js_text_template.substitute({'python_data': json.dumps(network),
                                        'network': rand_str})
    network_plot = html_template.substitute({'js_text': js_text})
    d3_download = """ <script src="http://d3js.org/d3.v3.min.js" charset="utf-8"></script>"""

    network_plot = d3_download + network_plot

    HTML(network_plot)
    return HTML(network_plot)
