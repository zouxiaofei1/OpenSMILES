import React, { useCallback } from "react";
import { Editor } from "ketcher-react";
import { StandaloneStructServiceProvider } from "ketcher-standalone";

const structServiceProvider = new StandaloneStructServiceProvider();

export function KetcherApp() {
  const onInit = useCallback((ketcher) => {
    window.ketcher = ketcher;
    window.dispatchEvent(new Event("ketcher-ready"));
  }, []);

  return (
    <div style={{ width: "100%", height: "100%" }}>
      <Editor
        staticResourcesUrl={import.meta.env.BASE_URL}
        structServiceProvider={structServiceProvider}
        onInit={onInit}
      />
    </div>
  );
}
