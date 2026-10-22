"""Check an already-running local app via HTTP and a real Streamlit session."""

import asyncio
from pathlib import Path
import tempfile
from urllib.request import ProxyHandler, build_opener

import gdstk
from google.protobuf.json_format import MessageToDict
from streamlit.proto.BackMsg_pb2 import BackMsg
from streamlit.proto.ForwardMsg_pb2 import ForwardMsg
from websockets.asyncio.client import connect


async def receive_run(connection) -> list[dict]:
    messages=[]
    for _ in range(600):
        data=await connection.recv()
        if not isinstance(data,bytes):
            continue
        message=ForwardMsg()
        message.ParseFromString(data)
        record=MessageToDict(message)
        assert not record.get("delta",{}).get("newElement",{}).get("exception"),record
        messages.append(record)
        if message.WhichOneof("type")=="script_finished":
            return messages
    raise AssertionError("The app did not finish its render.")


def elements(messages):
    return [message.get("delta",{}).get("newElement",{}) for message in messages]


async def check_session() -> None:
    async with asyncio.timeout(40):
        async with connect(
            "ws://127.0.0.1:8501/_stcore/stream",
            subprotocols=["streamlit"], proxy=None, open_timeout=10,
        ) as connection:
            request = BackMsg()
            request.rerun_script.query_string = ""
            await connection.send(request.SerializeToString())
            initial=await receive_run(connection)
            assert "From circuit to qubit." in str(initial)
            assert "5.972 GHz" in str(initial)
            buttons={element["button"]["label"]:element["button"]["id"]
                     for element in elements(initial) if "button" in element}
            for label in ["Find a design","Generate chip layout"]:
                request=BackMsg()
                request.rerun_script.query_string=""
                trigger=request.rerun_script.widget_states.widgets.add()
                trigger.id=buttons[label]
                trigger.trigger_value=True
                await connection.send(request.SerializeToString())
                completed=await receive_run(connection)
                print("PASS: live form submitted:",label)
            downloads=[element["downloadButton"] for element in elements(completed)
                       if "downloadButton" in element]
            chip=next(element for element in downloads if element.get("label")=="Download chip GDS")
            url=chip["url"]
            if url.startswith("/"):
                url="http://127.0.0.1:8501"+url
            opener=build_opener(ProxyHandler({}))
            with opener.open(url,timeout=10) as response:
                data=response.read()
                assert response.status==200 and len(data)>1000
            with tempfile.TemporaryDirectory() as temp:
                path=Path(temp)/"downloaded.gds"
                path.write_bytes(data)
                library=gdstk.read_gds(path)
                assert {cell.name for cell in library.top_level()}=={"CHIP"}
                assert library.unit==1e-6
            print(f"PASS: downloaded and parsed {len(data)} chip GDS bytes from the running app.")


if __name__ == "__main__":
    opener = build_opener(ProxyHandler({}))
    with opener.open("http://127.0.0.1:8501/_stcore/health", timeout=10) as response:
        assert response.status == 200
        assert response.read().strip() == b"ok"
    print("PASS: HTTP health endpoint responds ok.")
    asyncio.run(check_session())
