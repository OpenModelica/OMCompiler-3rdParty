
package org.omc.graphstream;


import java.io.IOException;
import java.net.UnknownHostException;

import org.graphstream.graph.Graph;
import org.graphstream.graph.implementations.MultiGraph;
import org.graphstream.stream.netstream.packing.Base64Unpacker;
import org.graphstream.stream.netstream.NetStreamReceiver;
import org.graphstream.stream.thread.ThreadProxyPipe;

public class Receiver {

    public static void main(String[] args) throws UnknownHostException, IOException, InterruptedException {

        Graph g = new MultiGraph("G",false,true);
        g.display();
        NetStreamReceiver net = new NetStreamReceiver(2001);

        ThreadProxyPipe pipe = net.getDefaultStream();
        pipe.addSink(g);
        while (true) {
            pipe.pump();
            Thread.sleep(100);
        }
    }
}
